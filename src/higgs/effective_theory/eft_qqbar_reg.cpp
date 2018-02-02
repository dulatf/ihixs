
#include "constants.h"
#include<vector>
#include "as_series.h"
#include "higgs_eft.h"
#include "cppchaplin.h"
#include <complex>
using namespace std;

// quark - antiquark reg
namespace HEFT{
    
    double qqb_nlo_reg(const double& z, const double& L){
        return qqb_nlo_r_lz0(z,L);
    }
    
    
    double qqb_nlo_r_lz0(const double& z, const double& L)
    {
        const double zb=1.-z;
        return 1./z * 32./27. * pow(zb,3.);
    }
    
    
    double qqbar_nlo_r_L0(const double& z)
    {
        const double zb=1.-z;
        return 1./z * 32./27. * pow(zb,3.);
    }
    
    
    
    double qqb_nnlo_reg(const double& z, const double& L){
        return qqb_nnlo_r_lz0(z,L)
        + qqb_nnlo_r_lz1(z,L)
        + qqb_nnlo_r_lz2(z,L);
    }
    double qqb_nnlo_r_lz0(const double& z, const double& L)
    {
        return    qqb_nnlo_r_lz0_const(z,L)
        + qqb_nnlo_r_lz0_logz(z,L)
        + qqb_nnlo_r_lz0_logz_sq(z,L)
        + qqb_nnlo_r_lz0_logz_cube(z,L) ;
    }
    
    double qqb_nnlo_r_lz0_const(const double& z, const double& L)
    {
        const double nf = consts::nf;
        const double Nc = QCD::Nc;
        const double z2 = consts::z2;
        const double z3 = consts::z3;
        complex<double> res =-(pow(-1 + pow(Nc,2),2)*(2 + 2*z + pow(z,2))*z2*chaplin::HPL(-1,z))/(8.*pow(Nc,3)*z) - (pow(-1 + pow(Nc,2),2)*pow(2 + z,2)*z2*chaplin::HPL(1,z))/(2.*pow(Nc,2)*z) + (pow(-1 + pow(Nc,2),2)*(432*z + 216*pow(z,2) + 144*pow(Nc,2)*pow(z,3))*chaplin::HPL(-1,0,z))/(864.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),2)*(216 - 648*Nc + 864*L*Nc + 144*pow(Nc,2) - 864*z + 432*Nc*z + 864*L*Nc*z - 216*pow(Nc,2)*z + 864*pow(z,2) + 108*Nc*pow(z,2) + 216*L*Nc*pow(z,2) + 216*pow(Nc,2)*pow(z,2) - 360*pow(z,3))*chaplin::HPL(1,0,z))/(864.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),2)*(-432 - 432*z - 216*pow(z,2))*chaplin::HPL(-1,-1,0,z))/(864.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),2)*(648 + 648*z + 324*pow(z,2))*chaplin::HPL(-1,0,0,z))/(864.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),2)*(-2160*Nc - 2160*Nc*z - 540*Nc*pow(z,2))*chaplin::HPL(0,1,0,z))/(864.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),2)*(-2592*Nc - 2592*Nc*z - 648*Nc*pow(z,2))*chaplin::HPL(1,0,0,z))/(864.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),2)*(-1728*Nc - 1728*Nc*z - 432*Nc*pow(z,2))*chaplin::HPL(1,1,0,z))/(864.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),2)*(369 - 2835*Nc - 1377*L*Nc - 324*pow(L,2)*Nc + 1386*pow(Nc,2) + 396*L*pow(Nc,2) - 92*Nc*nf - 72*L*Nc*nf + 432*z + 144*L*z + 2538*Nc*z + 972*L*Nc*z + 216*pow(L,2)*Nc*z - 4698*pow(Nc,2)*z - 1332*L*pow(Nc,2)*z + 444*Nc*nf*z + 216*L*Nc*nf*z - 585*pow(z,2) - 144*L*pow(z,2) + 297*Nc*pow(z,2) + 405*L*Nc*pow(z,2) + 108*pow(L,2)*Nc*pow(z,2) + 4590*pow(Nc,2)*pow(z,2) + 1332*L*pow(Nc,2)*pow(z,2) - 516*Nc*nf*pow(z,2) - 216*L*Nc*nf*pow(z,2) - 216*pow(z,3) - 1278*pow(Nc,2)*pow(z,3) - 396*L*pow(Nc,2)*pow(z,3) + 164*Nc*nf*pow(z,3) + 72*L*Nc*nf*pow(z,3) - 288*z2 + 648*Nc*z2 + 864*L*Nc*z2 - 144*pow(Nc,2)*z2 + 864*z*z2 - 432*Nc*z*z2 + 864*L*Nc*z*z2 + 648*pow(Nc,2)*z*z2 - 540*pow(z,2)*z2 - 324*Nc*pow(z,2)*z2 + 216*L*Nc*pow(z,2)*z2 - 648*pow(Nc,2)*pow(z,2)*z2 + 144*pow(z,3)*z2 + 360*pow(Nc,2)*pow(z,3)*z2 - 432*z3 - 432*z*z3 - 216*pow(z,2)*z3))/(864.*pow(Nc,3)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res);
    }
    
    double qqb_nnlo_r_lz0_logz(const double& z, const double& L)
    {
        const double nf = consts::nf;
        const double Nc = QCD::Nc;
        const double z2 = consts::z2;
        complex<double> res =
        -(pow(-1 + pow(Nc,2),2)*(12*pow(Nc,2)*(-11 + (53 + 6*L)*z - 2*(25 + 3*L)*pow(z,2) + (22 + 4*L)*pow(z,3)) - 6*z*(32 + 25*z + 12*pow(z,2) + 4*L*(3 - 3*z + 2*pow(z,2))) + Nc*(18*pow(L,2)*pow(2 + z,2) - 18*L*(-12 + 8*z + 5*pow(z,2)) - 8*nf*(-3 + 15*z - 12*pow(z,2) + 4*pow(z,3)) - 3*(3*(-59 + 44*z + 29*pow(z,2)) + 12*pow(2 + z,2)*z2))))/(288.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),2)*pow(2 + z,2)*chaplin::HPL(1,0,z))/(2.*pow(Nc,2)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*log(z);
    }
    
    double qqb_nnlo_r_lz0_logz_sq(const double& z, const double& L)
    {
        const double Nc = QCD::Nc;
        complex<double> res =
        -(pow(-1 + pow(Nc,2),2)*(z*(15 + 18*z - 2*pow(z,2)) + pow(Nc,2)*z*(15 - 18*z + 10*pow(z,2)) + 3*Nc*(L*pow(2 + z,2) - 2*(-3 + 4*z + pow(z,2)))))/(48.*pow(Nc,3)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*pow(log(z),2.);
    }
    
    double qqb_nnlo_r_lz0_logz_cube(const double& z, const double& L)
    {
        const double Nc = QCD::Nc;
        complex<double> res =
        -(pow(-1 + pow(Nc,2),2)*pow(2 + z,2))/(48.*pow(Nc,2)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*pow(log(z),3.);
    }
    
    double qqb_nnlo_r_lz1(const double& z, const double& L)
    {
        const double nf = consts::nf;
        const double Nc = QCD::Nc;
        const double z2 = consts::z2;
        
        complex<double> res =
        -(pow(-1 + pow(Nc,2),2)*((1 - z)*(12*pow(Nc,2)*((11 + 2*L)*pow(-1 + z,2) - 2*z) - 12*(3 + 2*L*pow(-1 + z,2) - 8*z + 3*pow(z,2)) - Nc*(8*nf*pow(-1 + z,2) + 72*L*(3 + z) + 27*(17 + 5*z))) + 72*Nc*pow(2 + z,2)*z2))/(144.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),2)*(4 - 18*z + 18*pow(z,2) - 8*pow(z,3) + 2*pow(Nc,2)*(-2 + 9*z - 9*pow(z,2) + 4*pow(z,3)) + 3*Nc*(12 - 8*z - 5*pow(z,2) + 2*L*pow(2 + z,2)))*chaplin::HPL(0,z))/(24.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),2)*pow(2 + z,2)*pow(chaplin::HPL(0,z),2))/(8.*pow(Nc,2)*z) - (pow(-1 + pow(Nc,2),2)*pow(2 + z,2)*chaplin::HPL(1,0,z))/(2.*pow(Nc,2)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*log(1.-z);
    }
    
    double qqb_nnlo_r_lz2(const double& z, const double& L)
    {
        const double Nc = QCD::Nc;
        
        complex<double> res =
        -(pow(-1 + pow(Nc,2),2)*(-1 + z)*(-pow(-1 + z,2) + 3*pow(Nc,2)*pow(-1 + z,2) - 6*Nc*(3 + z)))/(12.*pow(Nc,3)*z) - (pow(-1 + pow(Nc,2),2)*pow(2 + z,2)*chaplin::HPL(0,z))/(4.*pow(Nc,2)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*pow(log(1.-z),2.);
    }
    
    // mathematica file for qqb_nnlo_lzbar*_lz*_L* : Higgs_falko_qqb_n2

    
    double qqb_nnlo_r_lz0_L0(const double& z)
    {
        const double Nc = QCD::Nc;
        complex<double> res =-(pow(-1 + Nc,2)*pow(1 + Nc,2)*(216*consts::z2 + 216*z*consts::z2 + 108*pow(z,2)*consts::z2)*chaplin::HPL(-1,z))/(864.*pow(Nc,3)*z) - (pow(-1 + Nc,2)*pow(1 + Nc,2)*(1728*Nc*consts::z2 + 1728*Nc*z*consts::z2 + 432*Nc*pow(z,2)*consts::z2)*chaplin::HPL(1,z))/(864.*pow(Nc,3)*z) - (pow(-1 + Nc,2)*pow(1 + Nc,2)*(-432*z - 216*pow(z,2) - 144*pow(Nc,2)*pow(z,3))*chaplin::HPL(-1,0,z))/(864.*pow(Nc,3)*z) - (pow(-1 + Nc,2)*pow(1 + Nc,2)*(432 + 432*z + 216*pow(z,2))*chaplin::HPL(-1,-1,0,z))/(864.*pow(Nc,3)*z) - (pow(-1 + Nc,2)*pow(1 + Nc,2)*(-648 - 648*z - 324*pow(z,2))*chaplin::HPL(-1,0,0,z))/(864.*pow(Nc,3)*z) - (pow(-1 + Nc,2)*pow(1 + Nc,2)*(2160*Nc + 2160*Nc*z + 540*Nc*pow(z,2))*chaplin::HPL(0,1,0,z))/(864.*pow(Nc,3)*z) - (pow(-1 + Nc,2)*pow(1 + Nc,2)*(2592*Nc + 2592*Nc*z + 648*Nc*pow(z,2))*chaplin::HPL(1,0,0,z))/(864.*pow(Nc,3)*z) - (pow(-1 + Nc,2)*pow(1 + Nc,2)*(1728*Nc + 1728*Nc*z + 432*Nc*pow(z,2))*chaplin::HPL(1,1,0,z))/(864.*pow(Nc,3)*z) - (pow(-1 + Nc,2)*pow(1 + Nc,2)*chaplin::HPL(1,0,z)*(-216 + 648*Nc - 144*pow(Nc,2) + 864*z - 432*Nc*z + 216*pow(Nc,2)*z - 864*pow(z,2) - 108*Nc*pow(z,2) - 216*pow(Nc,2)*pow(z,2) + 360*pow(z,3) - 1728*Nc*log(z) - 1728*Nc*z*log(z) - 432*Nc*pow(z,2)*log(z)))/(864.*pow(Nc,3)*z) - (pow(-1 + Nc,2)*pow(1 + Nc,2)*(-369 + 2835*Nc - 1386*pow(Nc,2) + 92*Nc*consts::nf - 432*z - 2538*Nc*z + 4698*pow(Nc,2)*z - 444*Nc*consts::nf*z + 585*pow(z,2) - 297*Nc*pow(z,2) - 4590*pow(Nc,2)*pow(z,2) + 516*Nc*consts::nf*pow(z,2) + 216*pow(z,3) + 1278*pow(Nc,2)*pow(z,3) - 164*Nc*consts::nf*pow(z,3) + 288*consts::z2 - 648*Nc*consts::z2 + 144*pow(Nc,2)*consts::z2 - 864*z*consts::z2 + 432*Nc*z*consts::z2 - 648*pow(Nc,2)*z*consts::z2 + 540*pow(z,2)*consts::z2 + 324*Nc*pow(z,2)*consts::z2 + 648*pow(Nc,2)*pow(z,2)*consts::z2 - 144*pow(z,3)*consts::z2 - 360*pow(Nc,2)*pow(z,3)*consts::z2 + 1593*Nc*log(z) - 396*pow(Nc,2)*log(z) + 72*Nc*consts::nf*log(z) - 576*z*log(z) - 1188*Nc*z*log(z) + 1908*pow(Nc,2)*z*log(z) - 360*Nc*consts::nf*z*log(z) - 450*pow(z,2)*log(z) - 783*Nc*pow(z,2)*log(z) - 1800*pow(Nc,2)*pow(z,2)*log(z) + 288*Nc*consts::nf*pow(z,2)*log(z) - 216*pow(z,3)*log(z) + 792*pow(Nc,2)*pow(z,3)*log(z) - 96*Nc*consts::nf*pow(z,3)*log(z) - 432*Nc*consts::z2*log(z) - 432*Nc*z*consts::z2*log(z) - 108*Nc*pow(z,2)*consts::z2*log(z) + 324*Nc*pow(log(z),2) + 270*z*pow(log(z),2) - 432*Nc*z*pow(log(z),2) + 270*pow(Nc,2)*z*pow(log(z),2) + 324*pow(z,2)*pow(log(z),2) - 108*Nc*pow(z,2)*pow(log(z),2) - 324*pow(Nc,2)*pow(z,2)*pow(log(z),2) - 36*pow(z,3)*pow(log(z),2) + 180*pow(Nc,2)*pow(z,3)*pow(log(z),2) + 72*Nc*pow(log(z),3) + 72*Nc*z*pow(log(z),3) + 18*Nc*pow(z,2)*pow(log(z),3) + 432*consts::z3 + 432*z*consts::z3 + 216*pow(z,2)*consts::z3))/(864.*pow(Nc,3)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res);
    }
    
    double qqb_nnlo_r_lz0_L1(const double& z)
    {
        const double Nc = QCD::Nc;
       
        complex<double> res =-(pow(-1 + Nc,2)*pow(1 + Nc,2)*(-864*Nc - 864*Nc*z - 216*Nc*pow(z,2))*chaplin::HPL(1,0,z))/(864.*pow(Nc,3)*z) - (pow(-1 + Nc,2)*pow(1 + Nc,2)*(1377*Nc - 396*pow(Nc,2) + 72*Nc*consts::nf - 144*z - 972*Nc*z + 1332*pow(Nc,2)*z - 216*Nc*consts::nf*z + 144*pow(z,2) - 405*Nc*pow(z,2) - 1332*pow(Nc,2)*pow(z,2) + 216*Nc*consts::nf*pow(z,2) + 396*pow(Nc,2)*pow(z,3) - 72*Nc*consts::nf*pow(z,3) - 864*Nc*consts::z2 - 864*Nc*z*consts::z2 - 216*Nc*pow(z,2)*consts::z2 + 648*Nc*log(z) - 216*z*log(z) - 432*Nc*z*log(z) + 216*pow(Nc,2)*z*log(z) + 216*pow(z,2)*log(z) - 270*Nc*pow(z,2)*log(z) - 216*pow(Nc,2)*pow(z,2)*log(z) - 144*pow(z,3)*log(z) + 144*pow(Nc,2)*pow(z,3)*log(z) + 216*Nc*pow(log(z),2) + 216*Nc*z*pow(log(z),2) + 54*Nc*pow(z,2)*pow(log(z),2)))/(864.*pow(Nc,3)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res);
    }

    double qqb_nnlo_r_lz0_L2(const double& z)
    {
        const double Nc = QCD::Nc;

        complex<double> res =-(pow(-1 + Nc,2)*pow(1 + Nc,2)*(324*Nc - 216*Nc*z - 108*Nc*pow(z,2) + 216*Nc*log(z) + 216*Nc*z*log(z) + 54*Nc*pow(z,2)*log(z)))/(864.*pow(Nc,3)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res);
    }
    
    double qqb_nnlo_r_lz1_L0(const double& z)
    {
        const double Nc = QCD::Nc;
        
        complex<double> res =-(pow(-1 + pow(Nc,2),2)*(-36 - 459*Nc + 132*pow(Nc,2) - 8*Nc*consts::nf + 132*z + 324*Nc*z - 420*pow(Nc,2)*z + 24*Nc*consts::nf*z - 132*pow(z,2) + 135*Nc*pow(z,2) + 420*pow(Nc,2)*pow(z,2) - 24*Nc*consts::nf*pow(z,2) + 36*pow(z,3) - 132*pow(Nc,2)*pow(z,3) + 8*Nc*consts::nf*pow(z,3) + 288*Nc*consts::z2 + 288*Nc*z*consts::z2 + 72*Nc*pow(z,2)*consts::z2))/(144.*pow(Nc,3)*z) - (pow(-1 + pow(Nc,2),2)*(-24 - 216*Nc + 24*pow(Nc,2) + 108*z + 144*Nc*z - 108*pow(Nc,2)*z - 108*pow(z,2) + 90*Nc*pow(z,2) + 108*pow(Nc,2)*pow(z,2) + 48*pow(z,3) - 48*pow(Nc,2)*pow(z,3))*chaplin::HPL(0,z))/(144.*pow(Nc,3)*z) - (pow(-1 + pow(Nc,2),2)*(-72*Nc - 72*Nc*z - 18*Nc*pow(z,2))*pow(chaplin::HPL(0,z),2))/(144.*pow(Nc,3)*z) - (pow(-1 + pow(Nc,2),2)*(288*Nc + 288*Nc*z + 72*Nc*pow(z,2))*chaplin::HPL(1,0,z))/(144.*pow(Nc,3)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*log(1.-z);
    }
    
    double qqb_nnlo_r_lz1_L1(const double& z)
    {
        const double Nc = QCD::Nc;
        
        complex<double> res =-(pow(-1 + pow(Nc,2),2)*(-24 - 216*Nc + 24*pow(Nc,2) + 72*z + 144*Nc*z - 72*pow(Nc,2)*z - 72*pow(z,2) + 72*Nc*pow(z,2) + 72*pow(Nc,2)*pow(z,2) + 24*pow(z,3) - 24*pow(Nc,2)*pow(z,3)))/(144.*pow(Nc,3)*z) - (pow(-1 + pow(Nc,2),2)*(-144*Nc - 144*Nc*z - 36*Nc*pow(z,2))*chaplin::HPL(0,z))/(144.*pow(Nc,3)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*log(1.-z);
    }
    
    double qqb_nnlo_r_lz2_L0(const double& z)
    {
        const double Nc = QCD::Nc;
        
        complex<double> res =-(pow(-1 + pow(Nc,2),2)*(-1 + z)*(-pow(-1 + z,2) + 3*pow(Nc,2)*pow(-1 + z,2) - 6*Nc*(3 + z)))/(12.*pow(Nc,3)*z) - (pow(-1 + pow(Nc,2),2)*pow(2 + z,2)*chaplin::HPL(0,z))/(4.*pow(Nc,2)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*pow(log(1.-z),2.);
    }
    

    //--------------------------- n3lo
    //: implementation of Bernhard's triple expansion (Nov 2017)
    //: we have here the full qqbar reg, including logs(mh^2/muf^2)
    //: i=0: gg   i=1: qg   i=2: qqbar   i=3: qq   i=4: qQ2
    double N3LORegEvaluator::n3lo_reg_complete_qqbar(const double& z, unsigned int log_muf_mh_squared_power) {
        
        //: main forking into three regions of [0,1]
        double res = 0.0;
        if (z<=1./13.) {
            const double x = sqrt(z);
            const double logx = log(z);
            for (int t = 0; t < NumZTerms; t++) {
                for (unsigned int p=0; p < 6; p++) {
                    //: [a_s][initial state][log_muf_mh_squared_power][x power][log(x) power]
                    //: [3: n3lo]  [0: gg] ....
                    res += pow(x,t-2) * pow(logx,p) * ZExp[3][2][log_muf_mh_squared_power][t][p];
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
                res += pow(x,t) * WExp[3][2][log_muf_mh_squared_power][t][0];
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
                    res += pow(x,t) * pow(logx,p) * ZbExp[3][2][log_muf_mh_squared_power][t][p];
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
    ///-------------
    double qqb_n3lo_reg(const double& z, const double& L){
        return qqbar_n3lo_r_lz0_full_series(z,L)
        + qqbar_n3lo_r_lz1_full_series(z,L)
        + qqbar_n3lo_r_lz2_full_series(z,L)
        + qqbar_n3lo_r_lz3_full_series(z,L)
        + qqbar_n3lo_r_lz4_full_series(z,L)
        + qqbar_n3lo_r_lz5_full_series(z,L)
        + LEqqbN3LOregFalko(z,L);
    }
    
    double qqb_n3lo_reg_no_Lf(const double& z, const double& L){
        return qqbar_n3lo_r_lz0_full_series(z,L)
        + qqbar_n3lo_r_lz1_full_series(z,L)
        + qqbar_n3lo_r_lz2_full_series(z,L)
        + qqbar_n3lo_r_lz3_full_series(z,L)
        + qqbar_n3lo_r_lz4_full_series(z,L)
        + qqbar_n3lo_r_lz5_full_series(z,L);
    }
    
    double qqb_n3lo_reg_no_Lf(const double& z,int truncation_order){
        return qqbar_n3lo_r_lzX_series(0,z,truncation_order)
        + qqbar_n3lo_r_lzX_series(1,z,truncation_order)
        + qqbar_n3lo_r_lzX_series(2,z,truncation_order)
        + qqbar_n3lo_r_lzX_series(3,z,truncation_order)
        + qqbar_n3lo_r_lzX_series(4,z,truncation_order)
        + qqbar_n3lo_r_lzX_series(5,z,truncation_order);
    }
    

    
    
    double LEqqbN3LOregFalko_L3(const double& z){
        const double nf = consts::nf;
        const double Nc = QCD::Nc;
        const double z2= consts::z2;

        
        complex<double> L3=(pow(-1 + pow(Nc,2),2)*(-8640*pow(Nc,2)*nf*(-3 + 2*z + pow(z,2)) - 1080*Nc*(3*(-5 + z + 4*pow(z,2)) + 6*pow(2 + z,2)*z2) + 360*pow(Nc,3)*(-1438 + 1425*z - 3*pow(z,2) + 16*pow(z,3) + 54*pow(2 + z,2)*z2)))/(207360.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(-1 + 3*pow(Nc,2))*pow(2 + z,2)*chaplin::HPL(1,0,z))/(32.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),2)*(144*pow(Nc,2)*nf*pow(2 + z,2) - 54*Nc*z*(12 + 5*z) + 18*pow(Nc,3)*(-496 - 284*z + 55*pow(z,2)))*log(z))/(6912.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*(Nc*(-36*z - 9*pow(z,2)) + 9*pow(Nc,3)*(16 - 12*z + 9*pow(z,2)))*pow(log(z),2))/(1152.*pow(Nc,4)*z) + log(1 - z)*(-(pow(-1 + pow(Nc,2),2)*(-216*Nc*(-3 + 2*z + pow(z,2)) + 648*pow(Nc,3)*(-3 + 2*z + pow(z,2))))/(3456.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(-18*Nc*pow(2 + z,2) + 54*pow(Nc,3)*pow(2 + z,2))*log(z))/(576.*pow(Nc,4)*z))
        ;
        check_imaginary_part(L3,__PRETTY_FUNCTION__);
        return real(L3);
    }
    
    double LEqqbN3LOregFalko_L2(const double& z){
        const double nf = consts::nf;
        const double Nc = QCD::Nc;
        const double z2= consts::z2;
        const double z3= consts::z3;
        
        complex<double> L2=(pow(-1 + pow(Nc,2),2)*(36*(30*(-20*z + 20*pow(z,2)) + 60*(24*z - 24*pow(z,2) + 16*pow(z,3))*z2) + 1440*pow(Nc,4)*(121 + pow(z,2)*(455 - 36*z2) + pow(z,3)*(-121 + 24*z2) + z*(-455 + 36*z2)) - 4*pow(Nc,2)*(1440*pow(nf,2)*pow(-1 + z,3) + 90*nf*(-615 + 388*z + 227*pow(z,2) + 48*pow(2 + z,2)*z2) + 3*(20*(-642*z + 642*pow(z,2)) + 60*(144*z - 144*pow(z,2) + 96*pow(z,3))*z2)) - 120*pow(Nc,3)*(-24*nf*(-22 + 73*z - 73*pow(z,2) + 22*pow(z,3)) + 18*(-602 - 362*z + 85*pow(z,2))*z2 + 2*(11737 - 12165*z + 80*pow(z,3) + pow(z,2)*(348 - 648*z3) - 972*z3)) + 120*Nc*(2032 - 88*pow(z,3) + 108*(-15 + 5*z + 4*pow(z,2))*z2 + 3*pow(z,2)*(-9 + 56*nf - 72*z3) - 216*z3 - 3*z*(639 + 56*nf + 288*z3))))/(207360.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*(-60*Nc*pow(2 + z,2) + 180*pow(Nc,3)*pow(2 + z,2))*z2*chaplin::HPL(1,z))/(192.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(-3 - 2*z + pow(z,2))*chaplin::HPL(-1,0,z))/(4.*Nc*z) + (pow(-1 + pow(Nc,2),2)*(pow(Nc,3)*(494 + 482*z - 49*pow(z,2)) - 6*Nc*(-3 + 7*z + 2*pow(z,2)) + 8*(1 + pow(z,3)) + 8*pow(Nc,4)*(1 + pow(z,3)) - 8*pow(Nc,2)*(nf*pow(2 + z,2) + 2*(1 + pow(z,3))))*chaplin::HPL(1,0,z))/(96.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*pow(-2 + z,2)*chaplin::HPL(0,-1,0,z))/(8.*Nc*z) + (pow(-1 + pow(Nc,2),2)*(72*Nc + 36*pow(Nc,3)*(2 - 12*z + 3*pow(z,2)))*chaplin::HPL(0,1,0,z))/(192.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(-9*Nc*(2 + 4*z + pow(z,2)) + 18*pow(Nc,3)*(7 + 3*pow(z,2)))*chaplin::HPL(1,0,0,z))/(48.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*(-60*Nc*pow(2 + z,2) + 180*pow(Nc,3)*pow(2 + z,2))*chaplin::HPL(1,1,0,z))/(192.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(-12*(30 - 66*z)*z - 4*pow(Nc,4)*(1476*z - 1584*pow(z,2) + 924*pow(z,3)) - 8*pow(Nc,2)*(z*(-783 + 891*z - 462*pow(z,2)) + 3*nf*(-188 - 56*z + pow(z,2))) - 6*Nc*(-286 + 6*(213 + 28*nf)*z - 21*(-51 + 8*nf)*pow(z,2) + 16*(-1 + 7*nf)*pow(z,3) + 216*(2 + 4*z + pow(z,2))*z2) + 6*pow(Nc,3)*(-9216 + 2*(6275 + 84*nf)*z + (655 - 168*nf)*pow(z,2) + 16*(15 + 7*nf)*pow(z,3) + 144*(9 - 8*z + 6*pow(z,2))*z2))*log(z))/(6912.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*(306*Nc*z + 12*z*(3 - 3*z + 4*pow(z,2)) + 12*pow(Nc,4)*z*(3 - 3*z + 4*pow(z,2)) - 6*pow(Nc,3)*(-328 - 201*z + 63*pow(z,2)) - 2*pow(Nc,2)*(18*nf*pow(2 + z,2) + 12*z*(3 - 3*z + 4*pow(z,2))))*pow(log(z),2))/(1152.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*(-30*Nc*z*(4 + z) + 6*pow(Nc,3)*(48 - 60*z + 37*pow(z,2)))*pow(log(z),3))/(2304.*pow(Nc,4)*z) + pow(log(1 - z),2)*(-(pow(-1 + pow(Nc,2),2)*(4*pow(-1 + z,3) - 8*pow(Nc,2)*pow(-1 + z,3) + 4*pow(Nc,4)*pow(-1 + z,3) + 15*Nc*(-3 + 2*z + pow(z,2)) - 45*pow(Nc,3)*(-3 + 2*z + pow(z,2))))/(48.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*(-60*Nc*pow(2 + z,2) + 180*pow(Nc,3)*pow(2 + z,2))*log(z))/(384.*pow(Nc,4)*z)) + log(1 - z)*(-(pow(-1 + pow(Nc,2),2)*(36*(-16*z + 16*pow(z,2)) - 24*pow(Nc,4)*(-77 + 255*z - 255*pow(z,2) + 77*pow(z,3)) - 24*pow(Nc,2)*(-1 + z)*(-77 + 202*z - 77*pow(z,2) + 24*nf*(3 + z)) - 3*Nc*(112*nf*pow(-1 + z,3) + 3*(-771 + 408*z + 363*pow(z,2) + 120*pow(2 + z,2)*z2)) + 3*pow(Nc,3)*(112*nf*pow(-1 + z,3) + 3*(-6289 + 6128*z + 97*pow(z,2) + 64*pow(z,3) + 360*pow(2 + z,2)*z2))))/(3456.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*(-60*Nc*pow(2 + z,2) + 180*pow(Nc,3)*pow(2 + z,2))*chaplin::HPL(1,0,z))/(192.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(pow(Nc,3)*(4908 + 1596*z - 942*pow(z,2)) + 48*z*(3 - 3*z + 2*pow(z,2)) + 36*Nc*(-15 + 5*z + 4*pow(z,2)) + 8*pow(Nc,4)*(18*z - 18*pow(z,2) + 12*pow(z,3)) - 4*pow(Nc,2)*(72*z - 72*pow(z,2) + 48*pow(z,3) + 12*nf*pow(2 + z,2)))*log(z))/(576.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(-9*Nc*(2 + 4*z + pow(z,2)) + 18*pow(Nc,3)*(7 + 3*pow(z,2)))*pow(log(z),2))/(96.*pow(Nc,4)*z))
        ;
        check_imaginary_part(L2,__PRETTY_FUNCTION__);
        return real(L2);
    }
    
    double LEqqbN3LOregFalko_L1(const double& z){
        const double nf = consts::nf;
        const double Nc = QCD::Nc;
        const double z2= consts::z2;
        const double z3= consts::z3;
        const double z4= consts::z4;
        
        complex<double> L1=(pow(-1 + pow(Nc,2),2)*(120*pow(Nc,4)*(11679 - 45947*z + 45071*pow(z,2) - 10803*pow(z,3) + 12*(-154 + 951*z - 984*pow(z,2) + 484*pow(z,3))*z2 + 1080*z*z3 - 1296*pow(z,2)*z3 + 1296*pow(z,3)*z3) + 3*Nc*(288065 - 343900*z + 41035*pow(z,2) + 14800*pow(z,3) - 343560*z2 + 304560*z*z2 + 303480*pow(z,2)*z2 - 3840*pow(z,3)*z2 - 67680*z3 + 73440*z*z3 - 68040*pow(z,2)*z3 + 80*nf*(24*(5 - 6*z + 2*pow(z,3))*z2 + 3*(-52 + 33*pow(z,3) + 96*z3 + 12*pow(z,2)*(27 + 4*z3) + z*(-305 + 96*z3))) + 172800*z4 + 293760*z*z4 + 73440*pow(z,2)*z4) - pow(Nc,3)*(8545465 - 8483940*z - 240165*pow(z,2) + 178640*pow(z,3) - 2153880*z2 + 5256720*z*z2 + 196920*pow(z,2)*z2 + 138240*pow(z,3)*z2 + 80*nf*(4504 - 17511*z + 17952*pow(z,2) - 4945*pow(z,3) + 72*(-7 + 36*z - 33*pow(z,2) + 16*pow(z,3))*z2) + 626400*z3 + 997920*z*z3 + 191160*pow(z,2)*z3 + 272160*z4 - 2501280*z*z4 + 521640*pow(z,2)*z4) + 36*(60*(6 + 38*z + 55*pow(z,2) + 42*pow(z,3))*z2 + 30*(-13 + pow(z,2)*(1 - 372*z3) - 8*z3 - z*(1 + 48*z3) + pow(z,3)*(13 + 136*z3)) + 2610*z*(2 + z)*z4) - 4*pow(Nc,2)*(160*pow(nf,2)*(-1 + z)*(23 - 88*z + 41*pow(z,2)) + 60*nf*(6*(176 + 200*z + 53*pow(z,2))*z2 + z*(2687 - 216*z3) + pow(z,2)*(1027 - 54*z3) - 6*(619 + 36*z3)) + 3*(60*(238 + 168*z - 351*pow(z,2) + 478*pow(z,3))*z2 + 20*(pow(z,2)*(6325 - 1530*z3) - 2*z*(2903 - 954*z3) + 15*pow(z,3)*(43 + 84*z3) + 12*(-97 + 129*z3)) + 7830*z*(2 + z)*z4))))/(207360.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*(144*pow(Nc,4)*pow(z,3)*z2 + 864*pow(Nc,3)*(-3 - 2*z + pow(z,2))*z2 - 48*Nc*nf*(2 + 2*z + pow(z,2))*z2 - 3*(6*z*(28 + 13*z)*z2 + 96*(2 + 2*z + pow(z,2))*z3) + pow(Nc,2)*(6*(88 + 172*z + 83*pow(z,2) - 24*pow(z,3))*z2 + 288*(2 + 2*z + pow(z,2))*z3))*chaplin::HPL(-1,z))/(576.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*(6*(3*pow(Nc,3)*(992 + 1116*z - 129*pow(z,2)) - 3*Nc*(36 + 252*z + 73*pow(z,2)) + 16*pow(Nc,4)*(7 - 6*z + 6*pow(z,2) + 4*pow(z,3)) + 8*(-4 + 24*z - 33*pow(z,2) + 16*pow(z,3)) - 4*pow(Nc,2)*(20 + 24*z - 42*pow(z,2) + 48*pow(z,3) + 9*nf*pow(2 + z,2)))*z2 + 36*(8*(2 + 2*z + pow(z,2)) - 8*pow(Nc,2)*(2 + 2*z + pow(z,2)) - Nc*(20 + 68*z + 17*pow(z,2)) + pow(Nc,3)*(148 + 4*z + 97*pow(z,2)))*z3)*chaplin::HPL(1,z))/(1152.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),3)*(2 + 2*z + pow(z,2))*z2*chaplin::HPL(-1,-1,z))/(2.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(4*pow(Nc,4)*z*(-3 + 6*z + 44*pow(z,2)) - 4*Nc*(22 + 3*(3 + 8*nf)*z + 3*(-3 + 4*nf)*pow(z,2) + 4*pow(z,3)) + pow(Nc,3)*(-1071 - 792*z + 279*pow(z,2) - 32*nf*pow(z,3)) + 6*(-9*pow(z,2) + 24*(2 + 2*z + pow(z,2))*z2) - 6*pow(Nc,2)*(-(z*(90 + 49*z)) + 24*(2 + 2*z + pow(z,2))*z2))*chaplin::HPL(-1,0,z))/(288.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),3)*(2 + 2*z + pow(z,2))*z2*chaplin::HPL(-1,1,z))/(8.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(12*pow(Nc,3)*pow(-2 + z,2) + z*(2 + z) - pow(Nc,2)*z*(2 + z))*z2*chaplin::HPL(0,-1,z))/(16.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*(24 + pow(Nc,2)*(4 - 132*z + 31*pow(z,2)))*z2*chaplin::HPL(0,1,z))/(16.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),3)*(2 + 2*z + pow(z,2))*z2*chaplin::HPL(1,-1,z))/(8.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(36*(4 + 4*z + 13*pow(z,2) + 4*pow(z,3)) + 8*pow(Nc,4)*(88 - 60*z + 42*pow(z,2) + 55*pow(z,3)) - 4*pow(Nc,2)*(-8 + 576*z - 459*pow(z,2) + 366*pow(z,3) + nf*(104 + 248*z + 77*pow(z,2))) + Nc*(1311 + 234*z + 747*pow(z,2) - 16*pow(z,3) + 16*nf*(-6 + 33*z - 30*pow(z,2) + 13*pow(z,3)) + 72*(16 + 28*z + 7*pow(z,2))*z2) - 2*pow(Nc,3)*(4408 - 1291*z + 770*pow(z,2) + 96*pow(z,3) + 8*nf*(4 + 3*z + 3*pow(z,3)) + 36*(48 - 20*z + 27*pow(z,2))*z2))*chaplin::HPL(1,0,z))/(576.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(-49 + 143*pow(Nc,2))*pow(2 + z,2)*z2*chaplin::HPL(1,1,z))/(32.*pow(Nc,3)*z) - (pow(-1 + pow(Nc,2),2)*(8*pow(Nc,4)*pow(z,3) - 3*z*(12 + 5*z) + 48*pow(Nc,3)*(-3 - 2*z + pow(z,2)) - 8*Nc*nf*(2 + 2*z + pow(z,2)) + pow(Nc,2)*(88 + 124*z + 59*pow(z,2) - 8*pow(z,3)))*chaplin::HPL(-1,-1,0,z))/(48.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(16*pow(Nc,4)*pow(z,3) - 3*z*(28 + 11*z) + 120*pow(Nc,3)*(-3 - 2*z + pow(z,2)) - 24*Nc*nf*(2 + 2*z + pow(z,2)) + pow(Nc,2)*(264 + 348*z + 165*pow(z,2) - 16*pow(z,3)))*chaplin::HPL(-1,0,0,z))/(96.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*(2*pow(Nc,4)*pow(z,3) - 3*z*(2 + z) + pow(Nc,2)*z*(6 + 3*z - 2*pow(z,2)) + 12*pow(Nc,3)*(-3 - 2*z + pow(z,2)))*chaplin::HPL(-1,1,0,z))/(12.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(4*pow(Nc,4)*pow(z,3) - 6*z*(2 + z) - 3*pow(Nc,3)*(15 + 6*z - 5*pow(z,2)) + Nc*(-4 + 6*z - 3*pow(z,2)) + 2*pow(Nc,2)*z*(6 + 3*z - 2*pow(z,2)))*chaplin::HPL(0,-1,0,z))/(24.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(Nc*(344 + 276*z - 147*pow(z,2)) + pow(Nc,3)*(-3436 - 1908*z - 63*pow(z,2)) + 8*pow(Nc,4)*(4 - 6*z + 3*pow(z,2) + 2*pow(z,3)) + 8*(4 - 21*z - 12*pow(z,2) + 4*pow(z,3)) + 8*pow(Nc,2)*(-8 + 27*z + 9*pow(z,2) - 6*pow(z,3) + 9*nf*pow(2 + z,2)))*chaplin::HPL(0,1,0,z))/(192.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*(2*pow(Nc,4)*pow(z,3) - 3*z*(2 + z) + pow(Nc,2)*z*(6 + 3*z - 2*pow(z,2)) + 12*pow(Nc,3)*(-3 - 2*z + pow(z,2)))*chaplin::HPL(1,-1,0,z))/(12.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(-4 + 6*z - 72*pow(z,2) + 28*pow(z,3) + pow(Nc,3)*(-370 - 204*z - 237*pow(z,2)) - 6*Nc*(-6 - 2*z + 7*pow(z,2)) + pow(Nc,4)*(20 - 6*z + 20*pow(z,3)) + 8*pow(Nc,2)*(-2 + 9*pow(z,2) - 6*pow(z,3) + 3*nf*pow(2 + z,2)))*chaplin::HPL(1,0,0,z))/(48.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*(3*pow(Nc,3)*(1088 + 1180*z - 161*pow(z,2)) - 3*Nc*(36 + 252*z + 73*pow(z,2)) + 16*pow(Nc,4)*(7 - 6*z + 6*pow(z,2) + 3*pow(z,3)) + 16*(-2 + 15*z - 15*pow(z,2) + 8*pow(z,3)) - 4*pow(Nc,2)*(9*nf*pow(2 + z,2) + 4*(5 + 9*z - 9*pow(z,2) + 11*pow(z,3))))*chaplin::HPL(1,1,0,z))/(192.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),3)*(2 + 2*z + pow(z,2))*chaplin::HPL(-1,-1,-1,0,z))/(2.*pow(Nc,4)*z) - (5*pow(-1 + pow(Nc,2),3)*(2 + 2*z + pow(z,2))*chaplin::HPL(-1,-1,0,0,z))/(8.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),3)*(2 + 2*z + pow(z,2))*chaplin::HPL(-1,-1,1,0,z))/(4.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),3)*(2 + 2*z + pow(z,2))*chaplin::HPL(-1,0,-1,0,z))/(4.*pow(Nc,4)*z) + (3*pow(-1 + pow(Nc,2),3)*(2 + 2*z + pow(z,2))*chaplin::HPL(-1,0,0,0,z))/(8.*pow(Nc,4)*z) - (3*pow(-1 + pow(Nc,2),3)*(2 + 2*z + pow(z,2))*chaplin::HPL(-1,0,1,0,z))/(8.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),3)*(2 + 2*z + pow(z,2))*chaplin::HPL(-1,1,-1,0,z))/(4.*pow(Nc,4)*z) - (3*pow(-1 + pow(Nc,2),3)*(2 + 2*z + pow(z,2))*chaplin::HPL(-1,1,0,0,z))/(8.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(4*pow(Nc,3)*pow(-2 + z,2) + z*(2 + z) - pow(Nc,2)*z*(2 + z))*chaplin::HPL(0,-1,-1,0,z))/(8.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*(10*pow(Nc,3)*pow(-2 + z,2) + 3*z*(2 + z) - 3*pow(Nc,2)*z*(2 + z))*chaplin::HPL(0,-1,0,0,z))/(16.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*pow(-2 + z,2)*chaplin::HPL(0,-1,1,0,z))/(2.*Nc*z) - (pow(-1 + pow(Nc,2),2)*pow(-2 + z,2)*chaplin::HPL(0,0,-1,0,z))/(4.*Nc*z) - (pow(-1 + pow(Nc,2),2)*(-24 - 20*z - 5*pow(z,2) + pow(Nc,2)*(104 + 52*z + 21*pow(z,2)))*chaplin::HPL(0,0,1,0,z))/(32.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),2)*pow(-2 + z,2)*chaplin::HPL(0,1,-1,0,z))/(2.*Nc*z) + (pow(-1 + pow(Nc,2),2)*(16 - 4*z - pow(z,2) + 4*pow(Nc,2)*(-1 - 28*z + 6*pow(z,2)))*chaplin::HPL(0,1,0,0,z))/(16.*pow(Nc,3)*z) - (pow(-1 + pow(Nc,2),2)*(24 + pow(Nc,2)*(20 - 148*z + 35*pow(z,2)))*chaplin::HPL(0,1,1,0,z))/(16.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),3)*(2 + 2*z + pow(z,2))*chaplin::HPL(1,-1,-1,0,z))/(4.*pow(Nc,4)*z) - (3*pow(-1 + pow(Nc,2),3)*(2 + 2*z + pow(z,2))*chaplin::HPL(1,-1,0,0,z))/(8.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*pow(-2 + z,2)*chaplin::HPL(1,0,-1,0,z))/(2.*Nc*z) + (pow(-1 + pow(Nc,2),2)*(-2*(3 + 8*z + 2*pow(z,2)) + pow(Nc,2)*(42 - 12*z + 23*pow(z,2)))*chaplin::HPL(1,0,0,0,z))/(8.*pow(Nc,3)*z) - (3*pow(-1 + pow(Nc,2),2)*(2 + pow(Nc,2)*(2 - 12*z + 3*pow(z,2)))*chaplin::HPL(1,0,1,0,z))/(4.*pow(Nc,3)*z) - (pow(-1 + pow(Nc,2),2)*(-8*(5 + 8*z + 2*pow(z,2)) + pow(Nc,2)*(212 + 44*z + 83*pow(z,2)))*chaplin::HPL(1,1,0,0,z))/(16.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),2)*(-49 + 143*pow(Nc,2))*pow(2 + z,2)*chaplin::HPL(1,1,1,0,z))/(32.*pow(Nc,3)*z) + ((pow(-1 + pow(Nc,2),2)*(-4*pow(Nc,4)*(-2904 - 3*pow(z,2)*(6667 - 720*z2) + 3*z*(7057 - 648*z2) + 4*pow(z,3)*(2270 - 396*z2)) - Nc*(-12492 + 68058*z + 5712*nf*z + 19545*pow(z,2) + 192*nf*pow(z,2) + 2464*pow(z,3) + 2784*nf*pow(z,3) + 36*(392 - 324*z + 57*pow(z,2))*z2 + 3456*z*z3 + 864*pow(z,2)*z3) + pow(Nc,3)*(-193320 - 4224*nf + 299966*z + 23376*nf*z - 16369*pow(z,2) - 20256*nf*pow(z,2) + 864*pow(z,3) + 8672*nf*pow(z,3) + 23040*z2 + 11952*z*z2 - 9756*pow(z,2)*z2 + 5184*z3 + 17280*z*z3 + 5616*pow(z,2)*z3) - 12*z*(48*(15 + 6*z + pow(z,2))*z2 - 3*(57 + 241*z + 58*pow(z,2) + 48*z3 + 24*z*z3)) - 8*pow(Nc,2)*(16*pow(nf,2)*(-3 + 15*z - 12*pow(z,2) + 4*pow(z,3)) + nf*(-2213 + 460*z + 724*pow(z,2) + 36*pow(2 + z,2)*z2) + z*(-5445 + 2835*z - 2167*pow(z,2) + 36*(-3 - 42*z + 20*pow(z,2))*z2 + 216*z3 + 108*z*z3))))/(6912.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(11*Nc - 2*nf)*pow(2 + z,2)*chaplin::HPL(1,0,z))/(6.*pow(Nc,2)*z))*log(z) - (pow(-1 + pow(Nc,2),2)*(12*z*(20 + 19*z + 12*pow(z,2)) + 2*pow(Nc,4)*z*(1062 - 1257*z + 748*pow(z,2)) - 2*pow(Nc,2)*(2*nf*(188 + 8*z + 5*pow(z,2)) + 3*z*(-46 - 381*z + 156*pow(z,2))) + Nc*(-286 + (1125 - 96*nf)*z - 3*(-519 + 136*nf)*pow(z,2) + 112*(-1 + nf)*pow(z,3) + 18*(24 + 52*z + 13*pow(z,2))*z2) + pow(Nc,3)*(9176 - 18385*z - 384*nf*z - 559*pow(z,2) + 408*nf*pow(z,2) - 288*pow(z,3) - 240*nf*pow(z,3) - 18*(40 - 44*z + 45*pow(z,2))*z2))*pow(log(z),2))/(1152.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*(-3*Nc*z*(-228 + 49*z) + pow(Nc,3)*(2624 + 1860*z - 309*pow(z,2)) + 32*z*(-9 - 6*z + pow(z,2)) + 24*pow(Nc,4)*z*(6 - 7*z + 8*pow(z,2)) - 8*pow(Nc,2)*(6*nf*pow(2 + z,2) + z*(-18 - 45*z + 28*pow(z,2))))*pow(log(z),3))/(2304.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*(-11*z*(4 + z) + pow(Nc,2)*(96 - 132*z + 79*pow(z,2)))*pow(log(z),4))/(1536.*pow(Nc,3)*z) + pow(log(1 - z),3)*((pow(-1 + pow(Nc,2),2)*(8*pow(-1 + z,3) - 32*pow(Nc,2)*pow(-1 + z,3) + 24*pow(Nc,4)*pow(-1 + z,3) + 49*Nc*(-3 + 2*z + pow(z,2)) - 143*pow(Nc,3)*(-3 + 2*z + pow(z,2))))/(96.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(-49 + 143*pow(Nc,2))*pow(2 + z,2)*log(z))/(192.*pow(Nc,3)*z)) + pow(log(1 - z),2)*(-(pow(-1 + pow(Nc,2),2)*(288*(-1 + 4*z - 4*pow(z,2) + pow(z,3)) + 96*pow(Nc,4)*(-22 + 71*z - 71*pow(z,2) + 22*pow(z,3)) + 8*pow(Nc,2)*(-1 + z)*(27*nf*(3 + z) - 4*(53 - 130*z + 53*pow(z,2))) + Nc*(128*nf*pow(-1 + z,3) + 27*(-265 + 148*z + 117*pow(z,2)) + 882*pow(2 + z,2)*z2) - pow(Nc,3)*(256*nf*pow(-1 + z,3) + 3*(-12791 + 12492*z + 171*pow(z,2) + 128*pow(z,3) + 858*pow(2 + z,2)*z2))))/(1152.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(-49 + 143*pow(Nc,2))*pow(2 + z,2)*chaplin::HPL(1,0,z))/(64.*pow(Nc,3)*z) - (pow(-1 + pow(Nc,2),2)*(3*pow(Nc,3)*(2240 + 412*z - 545*pow(z,2)) - 3*Nc*(420 - 4*z - 55*pow(z,2)) + 16*(-4 + 21*z - 21*pow(z,2) + 10*pow(z,3)) + 16*pow(Nc,4)*(-4 + 27*z - 27*pow(z,2) + 14*pow(z,3)) - 4*pow(Nc,2)*(9*nf*pow(2 + z,2) + 32*(-1 + 6*z - 6*pow(z,2) + 3*pow(z,3))))*log(z))/(384.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*(-8*(5 + 8*z + 2*pow(z,2)) + pow(Nc,2)*(212 + 44*z + 83*pow(z,2)))*pow(log(z),2))/(64.*pow(Nc,3)*z)) + log(1 - z)*(-(pow(-1 + pow(Nc,2),2)*(-4*pow(Nc,4)*(-4603 + 16011*z - 15993*pow(z,2) + 4585*pow(z,3) - 72*(-4 + 27*z - 27*pow(z,2) + 15*pow(z,3))*z2) - 4*pow(Nc,2)*(2167 + 32*pow(nf,2)*pow(-1 + z,3) - 10548*z + 10593*pow(z,2) - 2212*pow(z,3) - 72*z2 + 1728*z*z2 - 2052*pow(z,2)*z2 + 1296*pow(z,3)*z2 + 3*nf*(-885 + 544*z + 341*pow(z,2) + 54*pow(2 + z,2)*z2) + 432*z3 + 432*z*z3 + 216*pow(z,2)*z3) - Nc*(-40448 + 31698*z + 8046*pow(z,2) + 704*pow(z,3) + 48*nf*(-23 + 95*z - 101*pow(z,2) + 29*pow(z,3)) + 20088*z2 + 1512*z*z2 - 2106*pow(z,2)*z2 + 2160*z3 + 7344*z*z3 + 1836*pow(z,2)*z3) + pow(Nc,3)*(-193858 + 201054*z - 5916*pow(z,2) - 1280*pow(z,3) + 16*nf*(-209 + 705*z - 723*pow(z,2) + 227*pow(z,3)) + 83376*z2 + 40392*z*z2 - 16902*pow(z,2)*z2 + 15984*z3 + 432*z*z3 + 10476*pow(z,2)*z3) + 36*(-36 + 29*pow(z,3) + 6*(4 - 4*z - 2*pow(z,2) + 4*pow(z,3))*z2 + 48*z3 + 8*pow(z,2)*(13 + 3*z3) + z*(-97 + 48*z3))))/(3456.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),3)*(2 + 2*z + pow(z,2))*z2*chaplin::HPL(-1,z))/(8.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(-49 + 143*pow(Nc,2))*pow(2 + z,2)*z2*chaplin::HPL(1,z))/(32.*pow(Nc,3)*z) - (pow(-1 + pow(Nc,2),2)*(2*pow(Nc,4)*pow(z,3) - 3*z*(2 + z) + pow(Nc,2)*z*(6 + 3*z - 2*pow(z,2)) + 12*pow(Nc,3)*(-3 - 2*z + pow(z,2)))*chaplin::HPL(-1,0,z))/(12.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*(3*pow(Nc,3)*(1112 + 1164*z - 169*pow(z,2)) - 3*Nc*(36 + 252*z + 73*pow(z,2)) + 16*pow(Nc,4)*(7 - 6*z + 6*pow(z,2) + 3*pow(z,3)) + 16*(-2 + 15*z - 15*pow(z,2) + 8*pow(z,3)) - 4*pow(Nc,2)*(9*nf*pow(2 + z,2) + 4*(5 + 9*z - 9*pow(z,2) + 11*pow(z,3))))*chaplin::HPL(1,0,z))/(192.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),3)*(2 + 2*z + pow(z,2))*chaplin::HPL(-1,-1,0,z))/(4.*pow(Nc,4)*z) - (3*pow(-1 + pow(Nc,2),3)*(2 + 2*z + pow(z,2))*chaplin::HPL(-1,0,0,z))/(8.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*pow(-2 + z,2)*chaplin::HPL(0,-1,0,z))/(2.*Nc*z) - (3*pow(-1 + pow(Nc,2),2)*(2 + pow(Nc,2)*(2 - 12*z + 3*pow(z,2)))*chaplin::HPL(0,1,0,z))/(4.*pow(Nc,3)*z) - (pow(-1 + pow(Nc,2),2)*(-8*(5 + 8*z + 2*pow(z,2)) + pow(Nc,2)*(212 + 44*z + 83*pow(z,2)))*chaplin::HPL(1,0,0,z))/(16.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),2)*(-49 + 143*pow(Nc,2))*pow(2 + z,2)*chaplin::HPL(1,1,0,z))/(32.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),2)*(12*z*(56 - 5*z + 24*pow(z,2)) + 8*pow(Nc,4)*(-77 + 477*z - 495*pow(z,2) + 220*pow(z,3)) - 4*pow(Nc,2)*(-154 + 1122*z - 1005*pow(z,2) + 512*pow(z,3) + nf*(284 + 128*z + 17*pow(z,2))) + pow(Nc,3)*(19705 - 25672*z - 1519*pow(z,2) - 480*pow(z,3) - 16*nf*(-7 + 36*z - 33*pow(z,2) + 14*pow(z,3)) - 3456*z2 + 1440*z*z2 - 1944*pow(z,2)*z2) + Nc*(-2867 + 2520*z + 2655*pow(z,2) - 32*pow(z,3) + 16*nf*(-7 + 36*z - 33*pow(z,2) + 14*pow(z,3)) + 1152*z2 + 2016*z*z2 + 504*pow(z,2)*z2))*log(z))/(576.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(-3*Nc*(30 - 32*z) + pow(Nc,3)*(730 + 236*z - 193*pow(z,2)) + 6*z*(-1 - 10*z + 4*pow(z,2)) + 2*pow(Nc,4)*z*(27 - 30*z + 20*pow(z,2)) - 8*pow(Nc,2)*(nf*pow(2 + z,2) + z*(6 - 15*z + 8*pow(z,2))))*pow(log(z),2))/(96.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(-2*(3 + 8*z + 2*pow(z,2)) + pow(Nc,2)*(42 - 12*z + 23*pow(z,2)))*pow(log(z),3))/(48.*pow(Nc,3)*z))
        ;
        check_imaginary_part(L1,__PRETTY_FUNCTION__);
        return real(L1);
    }
    
    
    double LEqqbN3LOregFalko(const double& z, const double& L)
    {
        if (abs(L)<1e-15) return 0.0;
        
        //const double nf = consts::nf;
        //const double Nc = QCD::Nc;
        //const double z2= consts::z2;
        //const double z3= consts::z3;
        //const double z4= consts::z4;
        
        double L3=LEqqbN3LOregFalko_L3(z);
        double L2 =LEqqbN3LOregFalko_L2(z) ;
        double L1 =LEqqbN3LOregFalko_L1(z);
        return L3*pow(L,3)+L2*pow(L,2)+L1*L;
    }
    
    //series expansions of the coeffs of log(1-z)^a
    
    
    
    double qqbar_n3lo_r_lz5_full_series(const double&z, const double& L){
        return qqbar_n3lo_r_lzX_series(5,z,30);
    }
    double qqbar_n3lo_r_lz4_full_series(const double&z, const double& L){
        return qqbar_n3lo_r_lzX_series(4,z,30);
    }
    double qqbar_n3lo_r_lz3_full_series(const double&z, const double& L){
        return qqbar_n3lo_r_lzX_series(3,z,30);
    }
    double qqbar_n3lo_r_lz2_full_series(const double&z, const double& L){
        return qqbar_n3lo_r_lzX_series(2,z,30);
    }
    double qqbar_n3lo_r_lz1_full_series(const double&z, const double& L){
        return qqbar_n3lo_r_lzX_series(1,z,30);
    }
    double qqbar_n3lo_r_lz0_full_series(const double&z, const double& L){
    return qqbar_n3lo_r_lzX_series(0,z,30);
    }
    
    class qqbarregn3locoeffs{
    public:
        static double give(int logpow,int coeff){
            double coeffs_lzbar_4[38]={0,7.456790123456790,11.18518518518519,29.27572016460905,34.86831275720165,39.59094650205761,43.69218106995885,47.31405055849500,50.55420340975897,53.48365667254556,56.15567313345091,58.61114139632658,60.88207293392479,62.99396104581290,64.96741557482298,66.81932169339577,68.56367795441870,70.21221341798684,71.77484958356418,73.26005133157983,74.67509717664516,76.02628997345200,77.31912306628509,78.55841267272445,79.74840437720891,80.89285955702641,81.99512609707073,83.05819668903021,84.08475723410927,85.07722729341628,86.03779409996010,86.96844132115237,87.87097351284227,88.74703701524986,89.59813789331152,90.42565740841274,91.23086541754231,92.01493202383114};
            
            double coeffs_lzbar_3[38]={0,-24.39506172839506,-39.66666666666667,-162.7475994513032,-153.9903978052126,-150.1102057613169,-145.9919478737997,-141.1304529548445,-135.6183077881358,-129.6223793893415,-123.2891118669690,-116.7318163439084,-110.0349160721513,-103.2605604118399,-96.45435828559468,-89.64976465901896,-82.87130445434514,-76.13691374069287,-69.45963896916869,-62.84887585998744,-56.31127876142643,-49.85143320658905,-43.47235718089708,-37.17587752307445,-30.96291455321997,-24.83369869710879,-18.78793631828294,-12.82493732519305,-6.943713805376547,-1.143056552564976,4.578405378945094,10.22215822512451,15.78976723481809,21.28284832470734,26.70304594658350,32.05201583824931,37.33141162620885,42.54287447485162};
            
            double coeffs_lzbar_2[38]={0,52.48989659765444,121.1422523038891,546.2618564904341,430.1066451938190,395.2026181144459,377.0324405298500,365.0568233084491,356.3053913771740,349.6483246237922,344.5442194816692,340.6802738750671,337.8484812124997,335.8965518350945,334.7058692836664,334.1803637575418,334.2402510065253,334.8181520825851,335.8564881122440,337.3056186264331,339.1224481284359,341.2693474332316,343.7132980643741,346.4252013129097,349.3793126783485,352.5527740400670,355.9252233712881,359.4784668298178,363.1962015762671,367.0637802148219,371.0680096416980,375.1969785203624,379.4399087089695,383.7870268307719,388.2294528626109,392.7591031628795,397.3686057998701,402.0512263974831};
            
            double coeffs_lzbar_1[38]={0,-13.56178692930848,-122.8388722197111,-747.6312247659640,-396.2995850133996,-305.8893443304928,-259.4270710482849,-228.0364980560464,-204.0698912691121,-184.6143729813627,-168.2530486314696,-154.1706002024059,-141.8419261018424,-130.9025761803966,-121.0865288938561,-112.1926717396459,-104.0650365506721,-96.58030274201162,-89.63946220712312,-83.16202732082040,-77.08187601020918,-71.34419529116624,-65.90318674811350,-60.72031490086829,-55.76295098007842,-51.00331002125280,-46.41760895057055,-41.98539341183551,-37.68899495351269,-33.51308997373632,-29.44433883790545,-25.47108869554984,-21.58312729780220,-17.77147793731082,-14.02822776257958,-10.34638334311517,-6.719748611442501,-3.142821277271304};

            
            
            
            double coeffs_lzbar_0[38]={0,-37.70751638827406,64.83686656787450,370.6425096527391,99.62093965349224,86.61281001454716,95.42583710432983,108.1223670459033,121.7312234626627,135.5057108662802,149.1865885116251,162.6590124436565,175.8643824142595,188.7715922384495,201.3652161860106,213.6395762214541,225.5952988675531,237.2371314861455,248.5724875924029,259.6104529258706,270.3610983406224,280.8350032212888,291.0429258511062,300.9955773712216,310.7034691653876,320.1768124196336,329.4254547620839,338.4588431984799,347.2860056068501,355.9155452255437,364.3556441237051,372.6140727621328,380.6982035607208,388.6150269744242,396.3711690052814,403.9729093878405,411.4261999111355,418.736682504977};
            
            
            switch (logpow){
                case 5: return 0.0;break;
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
    
    double qqbar_n3lo_r_lzX_series(int logzbpow,const double& z, int m){
        const int trunc_order = 37;
        if (m>trunc_order){
            cout<<"\nError: the qqbar_n3lo coefficients are truncated to O("<<trunc_order<<" and you asked for O("<<m<<"))"<<endl;
            exit(EXIT_FAILURE);
        }
        double res=0.;
        for (int i=0;i<m+1;i++){
            res += qqbarregn3locoeffs::give(logzbpow,i)*intpow(1.-z,i);
        }
        return res*intpow(log(1.-z),logzbpow);
    }

    
    
    
    
    
}
