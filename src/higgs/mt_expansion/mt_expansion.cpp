#include "mt_expansion.h"
#include "constants.h"
#include <iostream>





namespace MTEXP {

    double log_expanded(const double& z){
        return -1. +z - pow(-1 + z,2)/2. + pow(-1 + z,3)/3. - pow(-1 + z,4)/4. + pow(-1 + z,5)/5. - pow(-1 + z,6)/6. + pow(-1 + z,7)/7. - pow(-1 + z,8)/8. + pow(-1 + z,9)/9. - pow(-1 + z,10)/10. + pow(-1 + z,11)/11. - pow(-1 + z,12)/12. + pow(-1 + z,13)/13. - pow(-1 + z,14)/14. + pow(-1 + z,15)/15. ;
    }


    double gg_nlo_delta(const double& rho){
        const double a0 = 11./2.+6*consts::z2;
        const double a2 = 34./135.;
        const double a4 = 3553./113400.;
        const double a6 = 917641./190512000.;
        
        return a0+pow(rho,2.)*a2+pow(rho,4.)*a4+pow(rho,6.)*a6;
    }
    double gg_nlo_D0(const double& rho){return 0.0;}
    double gg_nlo_D1(const double& rho){return 12.0;}
    double gg_nlo_D0_L(const double& rho){return -6.0;}
    
    
    
    double gg_nlo_reg_he_limit(const double& z, const double& rho){
        const double tau=4./rho;
        //: data from Robert Harlander et al, http://arxiv.org/abs/0912.2104
        double Bgg1Table[23] = {-0.8821, 2.9212, 5.0234, 6.5538, 7.765, 8.7693, 9.6279, 10.3781,
            11.0444, 11.6437, 12.1883, 12.6875, 13.1482, 13.576, 13.9752,
            14.3495, 14.7018, 15.0345, 15.3497, 15.6491, 15.9343, 16.2065, 16.467};

        const double Bgg1 = MTEXP::linear_interpolate(Bgg1Table,tau);
        
        return Bgg1/z;
    }
    

    
    //: note that this is a very specific quick and dirty linear interpolator
    //: that ASSUMES a specific grid for tau values
    double linear_interpolate(double* my_table, const double&tau){
        for (int i=0;i<23;i++){
            double xgrid=double(i+2)/2.;
            if ((tau>xgrid) and (tau<=xgrid+0.5)){
                const double a=(my_table[i+1]-my_table[i])/0.5;
                const double b=my_table[i]-a*xgrid;
                return a*tau+b;
            }
        }
        if (tau>12.){
            
            return 0.0;
        }
        std::cout<<"[MTEXP]: interpolation attempt failed because x_tau=4*mt^2/mh^2 "
        <<"seems out of bounds : x_tau = "<<tau<<" not in [1.,12.]"<<std::endl;
        exit(1);
    }

    
    
    
    

    
    //: nnlo coefficients
    double gg_nnlo_delta(const double& rho,const double& lmH, const double& lmtop){
        

        const double nl = 5.;
        const double z2=consts::z2;
        const double z3=consts::z3;
        
        return 79.15972222222223 + (27*lmH)/2. + (19*lmtop)/8. - (1189*nl)/144. - (11*lmH*nl)/6. + (2*lmtop*nl)/3. - (47437199*rho)/1.24416e6 - (lmH*rho)/8. + (883*lmtop*rho)/1080. + (14563*nl*rho)/48600. + (1441*lmH*nl*rho)/25920. - (281*lmtop*nl*rho)/2880. - (998645169149*pow(rho,2))/1.170505728e11 + (1663*lmH*pow(rho,2))/201600. + (4039*lmtop*pow(rho,2))/51840. + (4565713*nl*pow(rho,2))/2.85768e8 + (80231*lmH*nl*pow(rho,2))/2.17728e7 - (193927*lmtop*nl*pow(rho,2))/2.17728e7 - (1712964005499545249*pow(rho,3))/3.93289924608e16 + (612511*lmH*pow(rho,3))/3.6288e8 + (88077779*lmtop*pow(rho,3))/7.62048e9 + (8432587511*nl*pow(rho,3))/4.8009024e12 + (5473619*lmH*nl*pow(rho,3))/1.306368e10 - (111726613*lmtop*nl*pow(rho,3))/9.144576e10 + (133*z2)/2. + (33*lmH*z2)/2. - 18*pow(lmH,2)*z2 - (5*nl*z2)/3. - lmH*nl*z2 + (89*rho*z2)/45. + (7*log(2.)*rho*z2)/45. - (7*nl*rho*z2)/90. + (9677*pow(rho,2)*z2)/37800. + (857*log(2.)*pow(rho,2)*z2)/37800. - (857*nl*pow(rho,2)*z2)/75600. + (646571*pow(rho,3)*z2)/1.5876e7 + (17881*log(2.)*pow(rho,3)*z2)/4.536e6 - (17881*nl*pow(rho,3)*z2)/9.072e6 - (9*pow(z2,2))/20. - (165*z3)/4. - (171*lmH*z3)/2. + (5*nl*z3)/6. + (1909181*rho*z3)/55296. + (267179777*pow(rho,2)*z3)/3.538944e7 + (5756378217151*pow(rho,3)*z3)/1.585446912e11;
    }
    
    double gg_nnlo_D0(const double& rho,const double& lmH){
        const double nl=5.;
        const double z2=consts::z2;
        const double z3=consts::z3;

        return -33.666666666666664 - (133*lmH)/2. - (33*pow(lmH,2))/4. + (14*nl)/9. + (5*lmH*nl)/3. + (pow(lmH,2)*nl)/2. - (68*lmH*rho)/45. - (3553*lmH*pow(rho,2))/18900. - (917641*lmH*pow(rho,3))/3.1752e7 + 33*z2 + 45*lmH*z2 - 2*nl*z2 + (351*z3)/2.;
    }
    
    double gg_nnlo_D1(const double& rho,const double& lmH){
        const double nl=5.;
        const double z2=consts::z2;

        return 133 + 33*lmH + 36*pow(lmH,2) - (10*nl)/3. - 2*lmH*nl + (136*rho)/45. + (3553*pow(rho,2))/9450. + (917641*pow(rho,3))/1.5876e7 - 90*z2;
    }
    
    double gg_nnlo_D2(const double& rho,const double& lmH){
        const double nl=5.;
        return -33 - 108*lmH + 2*nl;
    }
    
    double gg_nnlo_D3(const double& rho,const double& lmH){
        //const double nl=5.;
        return 72.;
    }

}


