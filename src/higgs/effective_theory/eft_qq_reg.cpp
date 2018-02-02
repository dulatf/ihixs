#include "constants.h"
#include<vector>
#include "as_series.h"
#include "higgs_eft.h"
#include "cppchaplin.h"
#include <complex>
using namespace std;

// quark - quark reg (same flavor)
namespace HEFT{
    
    double qq_nnlo_reg(const double& z, const double& L){
        return qq_nnlo_r_lz0(z,L)
        +qq_nnlo_r_lz1(z,L)
        +qq_nnlo_r_lz2(z,L);
    }
    
    double qq_nnlo_r_lz0(const double& z, const double& L)
    {
        return    qq_nnlo_r_lz0_const(z,L)
        + qq_nnlo_r_lz0_logz(z,L)
        + qq_nnlo_r_lz0_logz_sq(z,L)
        + qq_nnlo_r_lz0_logz_cube(z,L) ;
    }
    
    double qq_nnlo_r_lz0_const(const double& z, const double& L)
    {
        const double Nc = QCD::Nc;
        const double z2 = consts::z2;
        const double z3 = consts::z3;
        complex<double> res =
        -((pow(-1 + pow(Nc,2),2)*pow(2 + z,2)*z2*chaplin::HPL(1,z))/(pow(Nc,2)*z)) + (pow(-1 + pow(Nc,2),2)*(20 - 3*pow(z,2) + 2*L*pow(2 + z,2))*chaplin::HPL(0,1,z))/(8.*pow(Nc,2)*z) + (pow(-1 + pow(Nc,2),2)*(168*Nc + 192*L*Nc + 48*Nc*z + 192*L*Nc*z - 24*Nc*pow(z,2) + 48*L*Nc*pow(z,2))*chaplin::HPL(1,0,z))/(96.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),2)*(96*Nc + 96*Nc*z + 24*Nc*pow(z,2))*chaplin::HPL(0,0,1,z))/(96.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),2)*(-144*Nc - 144*Nc*z - 36*Nc*pow(z,2))*chaplin::HPL(0,1,0,z))/(96.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),2)*(192*Nc + 192*Nc*z + 48*Nc*pow(z,2))*chaplin::HPL(0,1,1,z))/(96.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),2)*(24 - 192*Nc - 24*z - 192*Nc*z + 12*pow(z,2) - 48*Nc*pow(z,2))*chaplin::HPL(1,0,0,z))/(96.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),2)*(-384*Nc - 384*Nc*z - 96*Nc*pow(z,2))*chaplin::HPL(1,1,0,z))/(96.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),2)*(-15 - 315*Nc - 153*L*Nc - 36*pow(L,2)*Nc - 48*z + 282*Nc*z + 108*L*Nc*z + 24*pow(L,2)*Nc*z + 63*pow(z,2) + 33*Nc*pow(z,2) + 45*L*Nc*pow(z,2) + 12*pow(L,2)*Nc*pow(z,2) + 72*Nc*z2 + 96*L*Nc*z2 - 48*Nc*z*z2 + 96*L*Nc*z*z2 - 36*Nc*pow(z,2)*z2 + 24*L*Nc*pow(z,2)*z2 - 24*z3 + 24*z*z3 - 12*pow(z,2)*z3))/(96.*pow(Nc,3)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res);
    }
    
    double qq_nnlo_r_lz0_logz(const double& z, const double& L)
    {
        const double Nc = QCD::Nc;
        const double z2 = consts::z2;
        complex<double> res =
        -(pow(-1 + pow(Nc,2),2)*(6*z*(4 + 9*z) + Nc*(6*pow(L,2)*pow(2 + z,2) - 6*L*(-12 + 8*z + 5*pow(z,2)) - 3*(-59 + 44*z + 29*pow(z,2)) - 12*pow(2 + z,2)*z2)))/(96.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),2)*pow(2 + z,2)*chaplin::HPL(1,0,z))/(2.*pow(Nc,2)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*log(z);
    }
    
    double qq_nnlo_r_lz0_logz_sq(const double& z, const double& L)
    {
        const double Nc = QCD::Nc;
        complex<double> res =
        -(pow(-1 + pow(Nc,2),2)*(-(z*(2 + 3*z)) + Nc*(L*pow(2 + z,2) - 2*(-3 + 4*z + pow(z,2)))))/(16.*pow(Nc,3)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*pow(log(z),2.);
    }
    
    double qq_nnlo_r_lz0_logz_cube(const double& z, const double& L)
    {
        const double Nc = QCD::Nc;
        complex<double> res =
        -(pow(-1 + pow(Nc,2),2)*pow(2 + z,2))/(48.*pow(Nc,2)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*pow(log(z),3.);
    }
    
    double qq_nnlo_r_lz1(const double& z, const double& L)
    {
        const double Nc = QCD::Nc;
        const double z2 = consts::z2;
        
        complex<double> res =
        (pow(-1 + pow(Nc,2),2)*(153 + 72*L - 108*z - 48*L*z - 45*pow(z,2) - 24*L*pow(z,2) - 192*z2 - 192*z*z2 - 48*pow(z,2)*z2))/(48.*pow(Nc,2)*z) + (pow(-1 + pow(Nc,2),2)*(L*pow(2 + z,2) - 2*(-4 + z + pow(z,2)))*chaplin::HPL(0,z))/(2.*pow(Nc,2)*z) + (pow(-1 + pow(Nc,2),2)*pow(2 + z,2)*pow(chaplin::HPL(0,z),2))/(4.*pow(Nc,2)*z) - (pow(-1 + pow(Nc,2),2)*pow(2 + z,2)*chaplin::HPL(1,0,z))/(pow(Nc,2)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*log(1.-z);
    }
    
    double qq_nnlo_r_lz2(const double& z, const double& L)
    {
        const double Nc = QCD::Nc;
        
        complex<double> res =
        -(pow(-1 + pow(Nc,2),2)*(3 - 2*z - pow(z,2)))/(2.*pow(Nc,2)*z) - (pow(-1 + pow(Nc,2),2)*pow(2 + z,2)*chaplin::HPL(0,z))/(2.*pow(Nc,2)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*pow(log(1.-z),2.);
    }
    
    // mathematica file for qq_n2lo_lzbar*_lz*_L* : Higgs_falko_qq_n2

    double qq_n2lo_lz0_L0(const double& z){
        const double Nc = QCD::Nc;
        //const double z2 = consts::z2;
        //const double z3 = consts::z3;
        complex<double> res =(pow(-1 + Nc,2)*pow(1 + Nc,2)*(-384*Nc*consts::z2 - 384*Nc*z*consts::z2 - 96*Nc*pow(z,2)*consts::z2)*chaplin::HPL(1,z))/(96.*pow(Nc,3)*z) + (pow(-1 + Nc,2)*pow(1 + Nc,2)*(240*Nc - 36*Nc*pow(z,2))*chaplin::HPL(0,1,z))/(96.*pow(Nc,3)*z) + (pow(-1 + Nc,2)*pow(1 + Nc,2)*(96*Nc + 96*Nc*z + 24*Nc*pow(z,2))*chaplin::HPL(0,0,1,z))/(96.*pow(Nc,3)*z) + (pow(-1 + Nc,2)*pow(1 + Nc,2)*(-144*Nc - 144*Nc*z - 36*Nc*pow(z,2))*chaplin::HPL(0,1,0,z))/(96.*pow(Nc,3)*z) + (pow(-1 + Nc,2)*pow(1 + Nc,2)*(192*Nc + 192*Nc*z + 48*Nc*pow(z,2))*chaplin::HPL(0,1,1,z))/(96.*pow(Nc,3)*z) + (pow(-1 + Nc,2)*pow(1 + Nc,2)*(24 - 192*Nc - 24*z - 192*Nc*z + 12*pow(z,2) - 48*Nc*pow(z,2))*chaplin::HPL(1,0,0,z))/(96.*pow(Nc,3)*z) + (pow(-1 + Nc,2)*pow(1 + Nc,2)*(-384*Nc - 384*Nc*z - 96*Nc*pow(z,2))*chaplin::HPL(1,1,0,z))/(96.*pow(Nc,3)*z) + (pow(-1 + Nc,2)*pow(1 + Nc,2)*chaplin::HPL(1,0,z)*(168*Nc + 48*Nc*z - 24*Nc*pow(z,2) + 192*Nc*log(z) + 192*Nc*z*log(z) + 48*Nc*pow(z,2)*log(z)))/(96.*pow(Nc,3)*z) + (pow(-1 + Nc,2)*pow(1 + Nc,2)*(-15 - 315*Nc - 48*z + 282*Nc*z + 63*pow(z,2) + 33*Nc*pow(z,2) + 72*Nc*consts::z2 - 48*Nc*z*consts::z2 - 36*Nc*pow(z,2)*consts::z2 - 177*Nc*log(z) - 24*z*log(z) + 132*Nc*z*log(z) - 54*pow(z,2)*log(z) + 87*Nc*pow(z,2)*log(z) + 48*Nc*consts::z2*log(z) + 48*Nc*z*consts::z2*log(z) + 12*Nc*pow(z,2)*consts::z2*log(z) - 36*Nc*pow(log(z),2) + 12*z*pow(log(z),2) + 48*Nc*z*pow(log(z),2) + 18*pow(z,2)*pow(log(z),2) + 12*Nc*pow(z,2)*pow(log(z),2) - 8*Nc*pow(log(z),3) - 8*Nc*z*pow(log(z),3) - 2*Nc*pow(z,2)*pow(log(z),3) - 24*consts::z3 + 24*z*consts::z3 - 12*pow(z,2)*consts::z3))/(96.*pow(Nc,3)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res);
    }
    
    double qq_n2lo_lz0_L1(const double& z){
        const double Nc = QCD::Nc;
         complex<double> res = (pow(-1 + Nc,2)*pow(1 + Nc,2)*(96*Nc + 96*Nc*z + 24*Nc*pow(z,2))*chaplin::HPL(0,1,z))/(96.*pow(Nc,3)*z) + (pow(-1 + Nc,2)*pow(1 + Nc,2)*(192*Nc + 192*Nc*z + 48*Nc*pow(z,2))*chaplin::HPL(1,0,z))/(96.*pow(Nc,3)*z) + (pow(-1 + Nc,2)*pow(1 + Nc,2)*(-153*Nc + 108*Nc*z + 45*Nc*pow(z,2) + 96*Nc*consts::z2 + 96*Nc*z*consts::z2 + 24*Nc*pow(z,2)*consts::z2 - 72*Nc*log(z) + 48*Nc*z*log(z) + 30*Nc*pow(z,2)*log(z) - 24*Nc*pow(log(z),2) - 24*Nc*z*pow(log(z),2) - 6*Nc*pow(z,2)*pow(log(z),2)))/(96.*pow(Nc,3)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res);
    }
    
    double qq_n2lo_lz0_L2(const double& z){
        const double Nc = QCD::Nc;
         complex<double> res = (pow(-1 + Nc,2)*pow(1 + Nc,2)*(-36*Nc + 24*Nc*z + 12*Nc*pow(z,2) - 24*Nc*log(z) - 24*Nc*z*log(z) - 6*Nc*pow(z,2)*log(z)))/(96.*pow(Nc,3)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res);
    }
    
    double qq_n2lo_lz1_L0(const double& z){
        const double Nc = QCD::Nc;
        complex<double> res =(pow(-1 + pow(Nc,2),2)*(153 - 108*z - 45*pow(z,2) - 192*consts::z2 - 192*z*consts::z2 - 48*pow(z,2)*consts::z2))/(48.*pow(Nc,2)*z) + (pow(-1 + pow(Nc,2),2)*(192 - 48*z - 48*pow(z,2))*chaplin::HPL(0,z))/(48.*pow(Nc,2)*z) + (pow(-1 + pow(Nc,2),2)*(48 + 48*z + 12*pow(z,2))*pow(chaplin::HPL(0,z),2))/(48.*pow(Nc,2)*z) + (pow(-1 + pow(Nc,2),2)*(-192 - 192*z - 48*pow(z,2))*chaplin::HPL(1,0,z))/(48.*pow(Nc,2)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*log(1.-z);
    }
    
    double qq_n2lo_lz1_L1(const double& z){
        const double Nc = QCD::Nc;
        complex<double> res =(pow(-1 + pow(Nc,2),2)*(72 - 48*z - 24*pow(z,2)))/(48.*pow(Nc,2)*z) + (pow(-1 + pow(Nc,2),2)*(96 + 96*z + 24*pow(z,2))*chaplin::HPL(0,z))/(48.*pow(Nc,2)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*log(1.-z);
    }
    
    double qq_n2lo_lz2_L0(const double& z){
        const double Nc = QCD::Nc;
        complex<double> res =-(pow(-1 + pow(Nc,2),2)*(3 - 2*z - pow(z,2)))/(2.*pow(Nc,2)*z) - (pow(-1 + pow(Nc,2),2)*pow(2 + z,2)*chaplin::HPL(0,z))/(2.*pow(Nc,2)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*pow(log(1.-z),2.);
    }
    
    
    //--------------------------- n3lo
    //: implementation of Bernhard's triple expansion (Nov 2017)
    //: we have here the full qqbar reg, including logs(mh^2/muf^2)
    //: i=0: gg   i=1: qg   i=2: qqbar   i=3: qq   i=4: qQ2
    double N3LORegEvaluator::n3lo_reg_complete_qq(const double& z, unsigned int log_muf_mh_squared_power) {
        
        //: main forking into three regions of [0,1]
        double res = 0.0;
        if (z<=1./13.) {
            const double x = sqrt(z);
            const double logx = log(z);
            for (int t = 0; t < NumZTerms; t++) {
                for (unsigned int p=0; p < 6; p++) {
                    //: [a_s][initial state][log_muf_mh_squared_power][x power][log(x) power]
                    //: [3: n3lo]  [0: gg] ....
                    res += pow(x,t-2) * pow(logx,p) * ZExp[3][3][log_muf_mh_squared_power][t][p];
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
                res += pow(x,t) * WExp[3][3][log_muf_mh_squared_power][t][0];
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
                    res += pow(x,t) * pow(logx,p) * ZbExp[3][3][log_muf_mh_squared_power][t][p];
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
    
    
    double qq_n3lo_reg(const double& z, const double& L){
        return qq_n3lo_r_lz0_full_series(z,L)
        +qq_n3lo_r_lz1_full_series(z,L)
        +qq_n3lo_r_lz2_full_series(z,L)
        +qq_n3lo_r_lz3_full_series(z,L)
        +qq_n3lo_r_lz4_full_series(z,L)
        +qq_n3lo_r_lz5_full_series(z,L)
        +LEqqN3LOregFalko(z,L);
    }
    
    double qq_n3lo_reg_no_Lf(const double& z, const double& L){
        return qq_n3lo_r_lz0_full_series(z,L)
        +qq_n3lo_r_lz1_full_series(z,L)
        +qq_n3lo_r_lz2_full_series(z,L)
        +qq_n3lo_r_lz3_full_series(z,L)
        +qq_n3lo_r_lz4_full_series(z,L)
        +qq_n3lo_r_lz5_full_series(z,L);
    }
    

    
    
    double LEqqN3LOregFalko_L1(const double& z)
    {
        const double nf = consts::nf;
        const double Nc = QCD::Nc;
        const double z2= consts::z2;
        const double z3= consts::z3;
        const double z4= consts::z4;
    
        complex<double> res= -(pow(-1 + pow(Nc,2),2)*(-72*(60*(3 + z - 14*pow(z,2) + 3*pow(z,3))*z2 + 15*(-13 - 8*z3 + pow(z,3)*(13 + 40*z3) + z*(-43 + 72*z3) + pow(z,2)*(43 + 126*z3)) + 540*(-2 + z)*z*z4) - 3*Nc*(288065 - 343900*z + 41035*pow(z,2) + 14800*pow(z,3) - 343560*z2 + 304560*z*z2 + 303480*pow(z,2)*z2 - 3840*pow(z,3)*z2 - 67680*z3 + 73440*z*z3 - 68040*pow(z,2)*z3 + 1440*nf*(5 - 8*z*(-2 + z3) + 8*z3 + pow(z,2)*(-21 + 4*z3)) + 172800*z4 + 293760*z*z4 + 73440*pow(z,2)*z4) + pow(Nc,3)*(360*(-5983 + 14602*z + 547*pow(z,2) + 384*pow(z,3))*z2 + 5*(1709093 + 35728*pow(z,3) + 125280*z3 + 81*pow(z,2)*(-593 + 472*z3) + 324*z*(-5237 + 616*z3)) + 3240*(84 - 772*z + 161*pow(z,2))*z4) + 24*pow(Nc,2)*(10*nf*(6*(176 + 200*z + 53*pow(z,2))*z2 + z*(2687 - 216*z3) + pow(z,2)*(1027 - 54*z3) - 6*(619 + 36*z3)) + 3*(60*(3 + z - 14*pow(z,2) + 3*pow(z,3))*z2 + 15*(97 + z*(309 - 104*z3) + 168*z3 + pow(z,3)*(13 + 40*z3) + pow(z,2)*(-419 + 214*z3)) + 540*(-2 + z)*z*z4))))/(207360.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*(1 + z)*(18*pow(Nc,3)*(-3 + z) + pow(1 + z,2) - pow(Nc,2)*pow(1 + z,2))*z2*chaplin::HPL(-1,z))/(12.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*(-48*(2 - 2*z + pow(z,2))*z3 - 12*pow(Nc,2)*(6*nf*pow(2 + z,2)*z2 - 4*(2 - 2*z + pow(z,2))*z3) - Nc*(6*(36 + 252*z + 73*pow(z,2))*z2 + 12*(20 + 68*z + 17*pow(z,2))*z3) + pow(Nc,3)*(6*(992 + 1116*z - 129*pow(z,2))*z2 + 12*(148 + 4*z + 97*pow(z,2))*z3))*chaplin::HPL(1,z))/(384.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(1 + z)*(9*pow(Nc,3)*(-119 + 31*z) + 12*(3 + 10*z + 3*pow(z,2)) - 12*pow(Nc,2)*(3 + 10*z + 3*pow(z,2)) - 4*Nc*(22 - 13*z + 4*pow(z,2)))*chaplin::HPL(-1,0,z))/(288.*pow(Nc,4)*z) + (3*pow(-1 + pow(Nc,2),2)*pow(-2 + z,2)*z2*chaplin::HPL(0,-1,z))/(4.*Nc*z) - (pow(-1 + pow(Nc,2),2)*(24 + pow(Nc,2)*(4 - 132*z + 31*pow(z,2)))*z2*chaplin::HPL(0,1,z))/(16.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),2)*(12*(-3*z*(4 + 9*z) + 6*(2 - 2*z + pow(z,2))*z2) - 4*pow(Nc,2)*(-9*z*(4 + 9*z) + nf*(104 + 248*z + 77*pow(z,2)) + 18*(2 - 2*z + pow(z,2))*z2) + Nc*(1311 + 234*z + 747*pow(z,2) - 16*pow(z,3) + 72*(16 + 28*z + 7*pow(z,2))*z2) + 2*pow(Nc,3)*(-4408 + 1291*z - 770*pow(z,2) - 96*pow(z,3) - 36*(48 - 20*z + 27*pow(z,2))*z2))*chaplin::HPL(1,0,z))/(576.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(-49 + 143*pow(Nc,2))*pow(2 + z,2)*z2*chaplin::HPL(1,1,z))/(32.*pow(Nc,3)*z) - (pow(-1 + pow(Nc,2),2)*(1 + z)*(6*pow(Nc,3)*(-3 + z) + pow(1 + z,2) - pow(Nc,2)*pow(1 + z,2))*chaplin::HPL(-1,-1,0,z))/(6.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(1 + z)*(15*pow(Nc,3)*(-3 + z) + pow(1 + z,2) - pow(Nc,2)*pow(1 + z,2))*chaplin::HPL(-1,0,0,z))/(12.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*(-3 - 2*z + pow(z,2))*chaplin::HPL(-1,1,0,z))/(Nc*z) - (pow(-1 + pow(Nc,2),2)*(3*pow(Nc,3)*(15 + 6*z - 5*pow(z,2)) - 2*z*(3 + 3*z + 2*pow(z,2)) + 2*pow(Nc,2)*z*(3 + 3*z + 2*pow(z,2)) + Nc*(4 - 6*z + 3*pow(z,2)))*chaplin::HPL(0,-1,0,z))/(24.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(24*z*(2 + 3*z) + Nc*(344 + 276*z - 147*pow(z,2)) + pow(Nc,3)*(-3436 - 1908*z - 63*pow(z,2)) + 24*pow(Nc,2)*(3*nf*pow(2 + z,2) - z*(2 + 3*z)))*chaplin::HPL(0,1,0,z))/(192.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*(-3 - 2*z + pow(z,2))*chaplin::HPL(1,-1,0,z))/(Nc*z) + (pow(-1 + pow(Nc,2),2)*(3*z*(12 + 11*z) + 2*pow(Nc,3)*(-370 - 204*z - 237*pow(z,2)) + pow(Nc,2)*(88 - 124*z + 11*pow(z,2) + 48*nf*pow(2 + z,2)) - 2*Nc*(4*nf*(2 - 2*z + pow(z,2)) + 6*(-6 - 2*z + 7*pow(z,2))))*chaplin::HPL(1,0,0,z))/(96.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*(-36 - 252*z - 73*pow(z,2) - 12*Nc*nf*pow(2 + z,2) + pow(Nc,2)*(1088 + 1180*z - 161*pow(z,2)))*chaplin::HPL(1,1,0,z))/(64.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),2)*pow(-2 + z,2)*chaplin::HPL(0,-1,-1,0,z))/(2.*Nc*z) - (5*pow(-1 + pow(Nc,2),2)*pow(-2 + z,2)*chaplin::HPL(0,-1,0,0,z))/(8.*Nc*z) + (pow(-1 + pow(Nc,2),2)*pow(-2 + z,2)*chaplin::HPL(0,-1,1,0,z))/(2.*Nc*z) - (pow(-1 + pow(Nc,2),2)*pow(-2 + z,2)*chaplin::HPL(0,0,-1,0,z))/(4.*Nc*z) - (pow(-1 + pow(Nc,2),2)*(-24 - 20*z - 5*pow(z,2) + pow(Nc,2)*(104 + 52*z + 21*pow(z,2)))*chaplin::HPL(0,0,1,0,z))/(32.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),2)*pow(-2 + z,2)*chaplin::HPL(0,1,-1,0,z))/(2.*Nc*z) + (pow(-1 + pow(Nc,2),2)*(-((-2 + z)*z) + pow(Nc,2)*(-2 + z)*z - Nc*(-16 + 4*z + pow(z,2)) + 4*pow(Nc,3)*(-1 - 28*z + 6*pow(z,2)))*chaplin::HPL(0,1,0,0,z))/(16.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*(24 + pow(Nc,2)*(20 - 148*z + 35*pow(z,2)))*chaplin::HPL(0,1,1,0,z))/(16.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),2)*pow(-2 + z,2)*chaplin::HPL(1,0,-1,0,z))/(2.*Nc*z) + (pow(-1 + pow(Nc,2),2)*(-2 + 2*z - pow(z,2) + pow(Nc,2)*(2 - 2*z + pow(z,2)) - 2*Nc*(3 + 8*z + 2*pow(z,2)) + pow(Nc,3)*(42 - 12*z + 23*pow(z,2)))*chaplin::HPL(1,0,0,0,z))/(8.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*(-2 + 12*Nc + 2*z - pow(z,2) + pow(Nc,2)*(2 - 2*z + pow(z,2)) + 6*pow(Nc,3)*(2 - 12*z + 3*pow(z,2)))*chaplin::HPL(1,0,1,0,z))/(8.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*(-2*(2 - 2*z + pow(z,2)) + 2*pow(Nc,2)*(2 - 2*z + pow(z,2)) - 8*Nc*(5 + 8*z + 2*pow(z,2)) + pow(Nc,3)*(212 + 44*z + 83*pow(z,2)))*chaplin::HPL(1,1,0,0,z))/(16.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(-49 + 143*pow(Nc,2))*pow(2 + z,2)*chaplin::HPL(1,1,1,0,z))/(32.*pow(Nc,3)*z) + ((pow(-1 + pow(Nc,2),2)*(-(Nc*(-12492 + 68058*z - 1152*nf*z + 19545*pow(z,2) - 2592*nf*pow(z,2) + 2464*pow(z,3) + 36*(392 - 324*z + 57*pow(z,2))*z2 + 3456*z*z3 + 864*pow(z,2)*z3)) + pow(Nc,3)*(-193320 + 299966*z - 16369*pow(z,2) + 864*pow(z,3) - 36*(-640 - 332*z + 271*pow(z,2))*z2 + 5184*z3 + 17280*z*z3 + 5616*pow(z,2)*z3) + 24*z*(12*(9 + 12*z + 2*pow(z,2))*z2 + 3*(-49 - 12*z3 + z*(-70 + 6*z3))) - 8*pow(Nc,2)*(nf*(-2213 + 460*z + 724*pow(z,2) + 36*pow(2 + z,2)*z2) + 3*z*(12*(9 + 12*z + 2*pow(z,2))*z2 + 3*(39 - 12*z3 + 2*z*(64 + 3*z3))))))/(6912.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(11*Nc - 2*nf)*pow(2 + z,2)*chaplin::HPL(1,0,z))/(6.*pow(Nc,2)*z))*log(z) - (pow(-1 + pow(Nc,2),2)*(12*z*(7 - 2*z + 6*pow(z,2)) - 4*pow(Nc,2)*(nf*(188 + 8*z + 5*pow(z,2)) + 3*z*(51 + 64*z + 6*pow(z,2))) + Nc*(-286 + (1125 + 96*nf)*z - 9*(-173 - 16*nf)*pow(z,2) - 112*pow(z,3) + 18*(24 + 52*z + 13*pow(z,2))*z2) + pow(Nc,3)*(9176 - 18385*z - 559*pow(z,2) - 288*pow(z,3) - 18*(40 - 44*z + 45*pow(z,2))*z2))*pow(log(z),2))/(1152.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*(-3*Nc*z*(-228 + 49*z) + pow(Nc,3)*(2624 + 1860*z - 309*pow(z,2)) + 16*z*(3 + 9*z + 2*pow(z,2)) - 16*pow(Nc,2)*(3*nf*pow(2 + z,2) + z*(3 + 9*z + 2*pow(z,2))))*pow(log(z),3))/(2304.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*(-11*z*(4 + z) + pow(Nc,2)*(96 - 132*z + 79*pow(z,2)))*pow(log(z),4))/(1536.*pow(Nc,3)*z) + pow(log(1 - z),3)*(-(pow(-1 + pow(Nc,2),2)*(-49 + 143*pow(Nc,2))*(-3 + 2*z + pow(z,2)))/(96.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),2)*(-49 + 143*pow(Nc,2))*pow(2 + z,2)*log(z))/(192.*pow(Nc,3)*z)) + pow(log(1 - z),2)*((pow(-1 + pow(Nc,2),2)*(2385 - 1332*z - 1053*pow(z,2) - 72*Nc*nf*(-3 + 2*z + pow(z,2)) - 1176*z2 - 1176*z*z2 - 294*pow(z,2)*z2 + pow(Nc,2)*(-12791 + 12492*z + 171*pow(z,2) + 128*pow(z,3) + 858*pow(2 + z,2)*z2)))/(384.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),2)*(-49 + 143*pow(Nc,2))*pow(2 + z,2)*chaplin::HPL(1,0,z))/(64.*pow(Nc,3)*z) - (pow(-1 + pow(Nc,2),2)*(-420 + 4*z + 55*pow(z,2) - 12*Nc*nf*pow(2 + z,2) + pow(Nc,2)*(2240 + 412*z - 545*pow(z,2)))*log(z))/(128.*pow(Nc,3)*z) - (pow(-1 + pow(Nc,2),2)*(-8*(5 + 8*z + 2*pow(z,2)) + pow(Nc,2)*(212 + 44*z + 83*pow(z,2)))*pow(log(z),2))/(64.*pow(Nc,3)*z)) + log(1 - z)*(-(pow(-1 + pow(Nc,2),2)*(-(Nc*(-40448 + 31698*z + 8046*pow(z,2) + 704*pow(z,3) - 54*(-372 - 28*z + 39*pow(z,2))*z2 + 2160*z3 + 7344*z*z3 + 1836*pow(z,2)*z3)) + 108*(5 - 8*z*(-2 + z3) + 8*z3 + pow(z,2)*(-21 + 4*z3)) - 12*pow(Nc,2)*(nf*(-885 + 544*z + 341*pow(z,2) + 54*pow(2 + z,2)*z2) + 9*(5 - 8*z*(-2 + z3) + 8*z3 + pow(z,2)*(-21 + 4*z3))) + pow(Nc,3)*(-54*(-1544 - 748*z + 313*pow(z,2))*z2 + 2*(-96929 - 640*pow(z,3) + 7992*z3 + 3*z*(33509 + 72*z3) + 6*pow(z,2)*(-493 + 873*z3)))))/(3456.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(-49 + 143*pow(Nc,2))*pow(2 + z,2)*z2*chaplin::HPL(1,z))/(32.*pow(Nc,3)*z) - (pow(-1 + pow(Nc,2),2)*(-3 - 2*z + pow(z,2))*chaplin::HPL(-1,0,z))/(Nc*z) - (pow(-1 + pow(Nc,2),2)*(-36 - 252*z - 73*pow(z,2) - 12*Nc*nf*pow(2 + z,2) + pow(Nc,2)*(1112 + 1164*z - 169*pow(z,2)))*chaplin::HPL(1,0,z))/(64.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),2)*pow(-2 + z,2)*chaplin::HPL(0,-1,0,z))/(2.*Nc*z) - (3*pow(-1 + pow(Nc,2),2)*(2 + pow(Nc,2)*(2 - 12*z + 3*pow(z,2)))*chaplin::HPL(0,1,0,z))/(4.*pow(Nc,3)*z) - (pow(-1 + pow(Nc,2),2)*(-2*(2 - 2*z + pow(z,2)) + 2*pow(Nc,2)*(2 - 2*z + pow(z,2)) - 8*Nc*(5 + 8*z + 2*pow(z,2)) + pow(Nc,3)*(212 + 44*z + 83*pow(z,2)))*chaplin::HPL(1,0,0,z))/(16.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(-49 + 143*pow(Nc,2))*pow(2 + z,2)*chaplin::HPL(1,1,0,z))/(32.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),2)*(-36*z*(4 + 9*z) - 4*pow(Nc,2)*(-9*z*(4 + 9*z) + nf*(284 + 128*z + 17*pow(z,2))) + Nc*(-2867 + 2520*z + 2655*pow(z,2) - 32*pow(z,3) + 72*(16 + 28*z + 7*pow(z,2))*z2) + pow(Nc,3)*(19705 - 25672*z - 1519*pow(z,2) - 480*pow(z,3) - 72*(48 - 20*z + 27*pow(z,2))*z2))*log(z))/(576.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(-3*Nc*(30 - 32*z) + 6*z*(2 + 3*z) + pow(Nc,3)*(730 + 236*z - 193*pow(z,2)) - 2*pow(Nc,2)*(4*nf*pow(2 + z,2) + 3*z*(2 + 3*z)))*pow(log(z),2))/(96.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(-2*(3 + 8*z + 2*pow(z,2)) + pow(Nc,2)*(42 - 12*z + 23*pow(z,2)))*pow(log(z),3))/(48.*pow(Nc,3)*z))
        
                        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res);
    }
    
    double LEqqN3LOregFalko_L2(const double& z)
    {
        const double nf = consts::nf;
        const double Nc = QCD::Nc;
        const double z2= consts::z2;
        const double z3= consts::z3;

        
        complex<double> res=-(pow(-1 + pow(Nc,2),2)*(360*pow(Nc,2)*nf*(-615 + 388*z + 227*pow(z,2) + 48*pow(2 + z,2)*z2) + 120*pow(Nc,3)*(18*(-602 - 362*z + 85*pow(z,2))*z2 + 2*(11737 - 12165*z + 80*pow(z,3) + pow(z,2)*(348 - 648*z3) - 972*z3)) - 120*Nc*(2032 - 88*pow(z,3) + 108*(-15 + 5*z + 4*pow(z,2))*z2 - 216*z3 - 27*pow(z,2)*(1 + 8*z3) - 27*z*(71 + 32*z3))))/(207360.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*(-120*Nc*pow(2 + z,2)*z2 + 360*pow(Nc,3)*pow(2 + z,2)*z2)*chaplin::HPL(1,z))/(384.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(-3 + z)*(1 + z)*chaplin::HPL(-1,0,z))/(4.*Nc*z) + (pow(-1 + pow(Nc,2),2)*(-48*pow(Nc,2)*nf*pow(2 + z,2) - 36*Nc*(-3 + 7*z + 2*pow(z,2)) - 6*pow(Nc,3)*(-494 - 482*z + 49*pow(z,2)))*chaplin::HPL(1,0,z))/(576.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*pow(-2 + z,2)*chaplin::HPL(0,-1,0,z))/(8.*Nc*z) + (pow(-1 + pow(Nc,2),2)*(72*Nc + 36*pow(Nc,3)*(2 - 12*z + 3*pow(z,2)))*chaplin::HPL(0,1,0,z))/(192.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(-18*Nc*(2 + 4*z + pow(z,2)) + 36*pow(Nc,3)*(7 + 3*pow(z,2)))*chaplin::HPL(1,0,0,z))/(96.*pow(Nc,4)*z) - (5*pow(-1 + pow(Nc,2),2)*(-1 + 3*pow(Nc,2))*pow(2 + z,2)*chaplin::HPL(1,1,0,z))/(16.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),2)*(-24*pow(Nc,2)*nf*(-188 - 56*z + pow(z,2)) - 6*Nc*(-286 + 1278*z + 1071*pow(z,2) - 16*pow(z,3) + 216*(2 + 4*z + pow(z,2))*z2) + 6*pow(Nc,3)*(-9216 + 12550*z + 655*pow(z,2) + 240*pow(z,3) + 144*(9 - 8*z + 6*pow(z,2))*z2))*log(z))/(6912.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*(306*Nc*z - 36*pow(Nc,2)*nf*pow(2 + z,2) - 6*pow(Nc,3)*(-328 - 201*z + 63*pow(z,2)))*pow(log(z),2))/(1152.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*(-30*Nc*z*(4 + z) + 6*pow(Nc,3)*(48 - 60*z + 37*pow(z,2)))*pow(log(z),3))/(2304.*pow(Nc,4)*z) + pow(log(1 - z),2)*((5*pow(-1 + pow(Nc,2),2)*(-1 + 3*pow(Nc,2))*(-3 + 2*z + pow(z,2)))/(16.*pow(Nc,3)*z) - (5*pow(-1 + pow(Nc,2),2)*(-1 + 3*pow(Nc,2))*pow(2 + z,2)*log(z))/(32.*pow(Nc,3)*z)) + log(1 - z)*(-(pow(-1 + pow(Nc,2),2)*(-576*pow(Nc,2)*nf*(-3 + 2*z + pow(z,2)) - 9*Nc*(-771 + 408*z + 363*pow(z,2) + 120*pow(2 + z,2)*z2) + 9*pow(Nc,3)*(-6289 + 6128*z + 97*pow(z,2) + 64*pow(z,3) + 360*pow(2 + z,2)*z2)))/(3456.*pow(Nc,4)*z) - (5*pow(-1 + pow(Nc,2),2)*(-1 + 3*pow(Nc,2))*pow(2 + z,2)*chaplin::HPL(1,0,z))/(16.*pow(Nc,3)*z) + (pow(-1 + pow(Nc,2),2)*(-48*pow(Nc,2)*nf*pow(2 + z,2) + pow(Nc,3)*(4908 + 1596*z - 942*pow(z,2)) + 36*Nc*(-15 + 5*z + 4*pow(z,2)))*log(z))/(576.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(-9*Nc*(2 + 4*z + pow(z,2)) + 18*pow(Nc,3)*(7 + 3*pow(z,2)))*pow(log(z),2))/(96.*pow(Nc,4)*z))
        
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res);
    }
    
    double LEqqN3LOregFalko_L3(const double& z)
    {
        const double nf = consts::nf;
        const double Nc = QCD::Nc;
        const double z2= consts::z2;

        
        complex<double> res=-(pow(-1 + pow(Nc,2),2)*(8640*pow(Nc,2)*nf*(-3 + 2*z + pow(z,2)) + 1080*Nc*(3*(-5 + z + 4*pow(z,2)) + 6*pow(2 + z,2)*z2) - 360*pow(Nc,3)*(-1438 + 1425*z - 3*pow(z,2) + 16*pow(z,3) + 54*pow(2 + z,2)*z2)))/(207360.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(-18*Nc*pow(2 + z,2) + 54*pow(Nc,3)*pow(2 + z,2))*chaplin::HPL(1,0,z))/(576.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(144*pow(Nc,2)*nf*pow(2 + z,2) - 54*Nc*z*(12 + 5*z) + 18*pow(Nc,3)*(-496 - 284*z + 55*pow(z,2)))*log(z))/(6912.*pow(Nc,4)*z) - (pow(-1 + pow(Nc,2),2)*(Nc*(-36*z - 9*pow(z,2)) + 9*pow(Nc,3)*(16 - 12*z + 9*pow(z,2)))*pow(log(z),2))/(1152.*pow(Nc,4)*z) + log(1 - z)*(-(pow(-1 + pow(Nc,2),2)*(-216*Nc*(-3 + 2*z + pow(z,2)) + 648*pow(Nc,3)*(-3 + 2*z + pow(z,2))))/(3456.*pow(Nc,4)*z) + (pow(-1 + pow(Nc,2),2)*(-18*Nc*pow(2 + z,2) + 54*pow(Nc,3)*pow(2 + z,2))*log(z))/(576.*pow(Nc,4)*z))
        
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res);
    }
    
    double LEqqN3LOregFalko(const double& z, const double& L)
    {
        if (abs(L)<1e-15) return 0.0;
        
        
        
        return   LEqqN3LOregFalko_L3(z)*pow(L,3)
                +LEqqN3LOregFalko_L2(z)*pow(L,2)
                +LEqqN3LOregFalko_L1(z)*L;
    }
    
    double qq_n3lo_reg_no_Lf(const double& z, int truncation_order){
        return qq_n3lo_r_lzX_series(0,z,truncation_order)
        +qq_n3lo_r_lzX_series(1,z,truncation_order)
        +qq_n3lo_r_lzX_series(2,z,truncation_order)
        +qq_n3lo_r_lzX_series(3,z,truncation_order)
        +qq_n3lo_r_lzX_series(4,z,truncation_order)
        +qq_n3lo_r_lzX_series(5,z,truncation_order);
    }
    
    
    
    
    
    
    //series expansions of the coeffs of log(1-z)^a
    
    
    
    double qq_n3lo_r_lz5_full_series(const double&z, const double& L){
        return qq_n3lo_r_lzX_series(5,z,30);
    }
    double qq_n3lo_r_lz4_full_series(const double&z, const double& L){
        return qq_n3lo_r_lzX_series(4,z,30);
    }
    double qq_n3lo_r_lz3_full_series(const double&z, const double& L){
        return qq_n3lo_r_lzX_series(3,z,30);
    }
    double qq_n3lo_r_lz2_full_series(const double&z, const double& L){
        return qq_n3lo_r_lzX_series(2,z,30);
    }
    double qq_n3lo_r_lz1_full_series(const double&z, const double& L){
        return qq_n3lo_r_lzX_series(1,z,30);
    }
    double qq_n3lo_r_lz0_full_series(const double&z, const double& L){
        return qq_n3lo_r_lzX_series(0,z,30);
    }
    
    class qqregn3locoeffs{
    public:
        static double give(int logpow,int coeff){
            double coeffs_lzbar_4[38]={0,7.456790123456790,11.18518518518519,18.64197530864198,24.23456790123457,28.95720164609053,33.05843621399177,36.68030570252792,39.92045855379189,42.84991181657848,45.52192827748383,47.97739654035950,50.24832807795771,52.36021618984582,54.33367071885590,56.18557683742869,57.92993309845162,59.57846856201976,61.14110472759710,62.62630647561275,64.04135232067809,65.39254511748492,66.68537821031801,67.92466781675738,69.11465952124183,70.25911470105933,71.36138124110365,72.42445183306313,73.45101237814219,74.44348243744921,75.40404924399302,76.33469646518530,77.23722865687519,78.11329215928278,78.96439303734444,79.79191255244566,80.59712056157524,81.38118716786406};
            
            double coeffs_lzbar_3[38]={0,-24.39506172839506,-39.66666666666667,-69.84773662551440,-81.46913580246914,-86.77577503429355,-88.44791495198903,-87.77407435402144,-85.52612671817698,-82.19732223342012,-78.11492674710488,-73.50246225035680,-68.51629897822552,-63.26805342180203,-57.83875422801532,-52.28805951459777,-46.66040652739783,-40.98920701615245,-35.29976677487384,-29.61135371899093,-23.93868631333805,-18.29302024012655,-12.68295203139561,-7.115020333038816,-1.594160512614378,3.875948331346309,9.292616197599394,14.65393486522295,19.95861587862040,25.20586312123184,30.39527201204689,35.52674934666765,40.60044925754544,45.61672184021996,50.57607179078040,55.47912499927462,60.32660149750805,65.11929350559348};
            
            double coeffs_lzbar_2[38]={0,52.48989659765444,115.8829930446298,206.8914081608028,237.1672668230354,253.8531162868722,264.5069030369813,271.8876173852453,277.4772355715482,282.1103646724569,286.2659369197624,290.2220918949855,294.1409336850088,298.1160829185782,302.2000366144192,306.4202900753142,310.7890410605265,315.3091377960698,319.9777786145654,324.7888385461644,329.7343426388407,334.8054010687596,339.9928005676627,345.2873742563688,350.6802276193454,356.1628707337101,361.7272894006198,367.3659766319043,373.0719386814969,378.8386850485087,384.6602087269864,390.5309608754629,396.4458226688779,402.4000761444816,408.3893752103330,414.4097175498466,420.4574178617595,426.5290826770403};
            
            double coeffs_lzbar_1[38]={0,-13.56178692930848,-100.4438104913160,-197.0289711735459,-201.4950538493328,-196.7023291346468,-189.7294769031933,-181.9018069878534,-174.0130498581842,-166.4410382184307,-159.3299303574603,-152.7088770637718,-146.5548856088234,-140.8240779993118,-135.4667277121877,-130.4341974626872,-125.6818753788213,-121.1701591628293,-116.8645149676797,-112.7351255776926,-108.7563839589466,-104.9063567710347,-101.1662753747808,-97.52007784628723,-93.95400856230057,-90.45627383480591,-87.01674856042556,-83.62672766328424,-80.27871606398868,-76.96625138083011,-73.68375424067451,-70.42640178473926,-67.19002062219491,-63.97099608066889,-60.76619511859561,-57.57290070247191,-54.38875581987755,-51.21171560569928};


            
            double coeffs_lzbar_0[38]={0,-39.01478263128528,16.21497910036776,49.52495968771541,45.64789732873069,49.19264772033760,56.53442969152410,65.48870254849564,75.30840402910120,85.59358020635574,96.08904207541665,106.6255137945448,117.0906777588666,127.4109532158581,137.5392482430913,147.4466110524096,157.1165165284698,166.5409428062287,175.7176657371448,184.6483846030208,193.3374179832537,201.7907930498192,210.0156081020876,218.0195860918832,225.8107624548055,233.3972678960460,240.7871786124447,247.9884145729458,255.0086721226187,261.8553811169570,268.5356795664871,275.0564007367166,281.4240690498371,287.6449021401005,293.7248171402318,299.6694398022335,305.4841154389733,311.1739209529166};
            
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
    
    double qq_n3lo_r_lzX_series(int logzbpow,const double& z, int m){
        const int trunc_order = 38;
        if (m>trunc_order){
            cout<<"\nError: the qg_n3lo coefficients are truncated to O("<<trunc_order<<" and you asked for O("<<m<<"))"<<endl;
            exit(EXIT_FAILURE);
        }
        double res=0.;
        for (int i=0;i<m+1;i++){
            res += qqregn3locoeffs::give(logzbpow,i)*intpow(1.-z,i);
        }
        return res*intpow(log(1.-z),logzbpow);
    }

    
    
    

}
