#include "cppchaplin.h"
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <iostream> /* cout */
#include "src/tools/constants.h"


namespace chaplin{


    complex<double> HPL2at0(int n1);
    complex<double> HPL2at1(int n1);
    complex<double> HPL2atm1(int n1);
    complex<double> HPL2ar1(int n1,const complex<double> & xx);
    complex<double> HPL2ar0(int n1,const complex<double> & xx);
    complex<double> HPL2arm1(int n1,const complex<double> & xx);
    complex<double> HPL2else(int n1,const complex<double> & xx);
    
    complex<double> HPL(int n1, int n2, const complex<double>& x){
        //chaplin::check_index(n1);
        //chaplin::check_index(n2);
        
        
        int region = chaplin::find_region(x);
                
        int index_in_base_3 = n2+1 + 3 * (n1+1);
        
        switch (region) {
            case 0:
                return chaplin::HPL2at0(index_in_base_3);
                break;
            case 1:
                return chaplin::HPL2at1(index_in_base_3);
                break;
            case 2:
                return chaplin::HPL2atm1(index_in_base_3);
                break;
            case 3:
                return chaplin::HPL2ar1(index_in_base_3,x);
                break;
            case 4:
                return chaplin::HPL2arm1(index_in_base_3,x);
                break;
            case 5:
                return chaplin::HPL2ar0(index_in_base_3,x);
                break;
            default:
                return chaplin::HPL2else(index_in_base_3,x);
                break;
        }
    }
    
    
    
    bool imaginary_part_should_be_zero_in_hpl2(int index, const complex<double>& x) {
        double xre = x.real();
        int n1 = index/3 -1;
        int n2 = index%3 -1;
        vector<int> vint = {n1,n2};
        if (index != chaplin::find_index_in_base_3(vint)) {
            cout << "shit in inverse index function " << endl;
        }
        if (n2 == 0 and xre > 0.0) {
            if (xre < 1.0) return true;
        }
        else if (n2 == 1 and xre < 1.0) {
            if (n1 != -1) return true;
            else if (xre > -1.0) return true;
        }
        else if (n2 == -1 and xre > -1.0) {
            if (n1 != 1) return true;
            else if (xre < 1.0) return true;
        }
        return false;
    }
    
    
    complex<double> HPL2at0(int index){
        if (index == 4) {
            cout << "chaplin error: call for HPL2(0,0,0.0)"
            << " which is divergent! "
            << endl;
            exit(EXIT_FAILURE);
        }
        else return complex<double>(0.0,0.0);
    }
    
    complex<double> HPL2at1(int index){
        complex<double> ll2_sq ( pow( log(2.0) , 2.0 )  ,  0.0);
        complex<double> pi_sq (consts::pi_square ,0.0 );
        const complex<double> zero(0.0,0.0);
        switch (index) {
            case 0: return ll2_sq/2.0;
                break;
            case 1: return -pi_sq/12.0;
                break;
            case 2: return pi_sq/12.0 - ll2_sq/2.0;
                break;
            case 3: return pi_sq/12.0;
                break;
            case 4: return zero;
                break;
            case 5: return pi_sq/6.0;
                break;
            case 6: {cout << "chaplin error: call for HPL2(1,0,1.0)"
                << " which is divergent! "
                << endl;
                exit(EXIT_FAILURE);
            }
                break;
            case 7: return -pi_sq/6.0;
                break;
            case 8: {cout << "chaplin error: call for HPL2(1,1,1.0)"
                << " which is divergent! "
                << endl;
                exit(EXIT_FAILURE);
            }
                break;
            default: {cout << "chaplin error: index in HPL2at1 out of bounds"
                << " : should be 0-8 and it is " << index << endl;
                exit(EXIT_FAILURE);
            }
        }
    }
    
    complex<double> HPL2atm1(int index){
        complex<double> ll2_sq ( pow( log(2.0) , 2.0 )  ,  0.0);
        complex<double> pi_sq (consts::pi_square ,0.0 );
        switch (index) {
            case 0: {cout << "chaplin error: call for HPL2(-1,-1,-1.0)"
                << " which is divergent! "
                << endl;
                exit(EXIT_FAILURE);}
                break;
            case 1: {cout << "chaplin error: call for HPL2(-1,0,-1.0)"
                << " which is divergent! "
                << endl;
                exit(EXIT_FAILURE);}
                break;
            case 2: {cout << "chaplin error: call for HPL2(-1,1,-1.0)"
                << " which is divergent! "
                << endl;
                exit(EXIT_FAILURE);}
                break;
            case 3:return - pi_sq / 6.0;
                break;
            case 4:return - pi_sq / 2.0;
                break;
            case 5:return  - pi_sq / 12.0;
                break;
            case 6:return pi_sq / 12.0 - ll2_sq / 2.0;
                break;
            case 7:return pi_sq / 12.0 - complex<double>(0.0,consts::Pi * log(2.));
                break;
            case 8:return ll2_sq / 2.0;
                break;
            default: {cout << "chaplin error: index in HPL2atm1 out of bounds"
                << " : should be 0-8 and it is " << index << endl;
                exit(EXIT_FAILURE);
            }
        }
    }
    
    
    complex<double> HPL2ar1(int index,const complex<double> & xx){
        
        const complex<double> ll2 ( log(2.0) ,  0.0);
        const complex<double> ll2_sq ( pow( log(2.0) , 2.0 )  ,  0.0);
        const complex<double> pi_sq (consts::pi_square ,0.0 );
        const double pi=consts::Pi;
        complex<double> x=xx;
        
        // adding i*epsilon
        bool imaginary_part_was_modified = false;
        if ( x.imag() == 0.0 ) {
            x = x + complex<double> (0.0, 1e-60);
            imaginary_part_was_modified = true;
        }
        complex<double> zp = 1.0-x;
        complex<double> llzp = log(zp);
        complex<double> res(0.0,0.0);
        
        
        switch (index) {
            case 0:  res = 0.5*pow(ll2,2.0) - 0.5*ll2*zp + (0.125 - 0.125*ll2)*pow(zp,2.0) + (0.0625 - 0.041666666666666664*ll2)*pow(zp,3.0) + (0.028645833333333332 - 0.015625*ll2)*pow(zp,4.0) + (0.013020833333333332 - 0.00625*ll2)*pow(zp,5.0) + (0.005946180555555556 - 0.0026041666666666665*ll2)*pow(zp,6.0);
                break;
            case 1:  res = -0.08333333333333333*pow(pi,2.0) + 0.25*pow(zp,2.0) + 0.16666666666666666*pow(zp,3.0) + 0.10416666666666666*pow(zp,4.0) + 0.06666666666666667*pow(zp,5.0) + 0.044444444444444446*pow(zp,6.0);
                break;
            case 2:  res = -0.5*pow(ll2,2.0) + 0.08333333333333333*pow(pi,2.0) + (-0.5 + 0.5*llzp)*zp + (-0.0625 + 0.125*llzp)*pow(zp,2.0) + (-0.013888888888888888 + 0.041666666666666664*llzp)*pow(zp,3.0) + (-0.00390625 + 0.015625*llzp)*pow(zp,4.0) + (-0.00125 + 0.00625*llzp)*pow(zp,5.0) + (-0.00043402777777777775 + 0.0026041666666666665*llzp)*pow(zp,6.0);
                break;
            case 3:  res = 0.08333333333333333*pow(pi,2.0) - 1.*ll2*zp + (0.25 - 0.5*ll2)*pow(zp,2.0) + (0.20833333333333331 - 0.3333333333333333*ll2)*pow(zp,3.0) + (0.16666666666666666 - 0.25*ll2)*pow(zp,4.0) + (0.13645833333333332 - 0.2*ll2)*pow(zp,5.0) + (0.11475694444444445 - 0.16666666666666666*ll2)*pow(zp,6.0);
                break;
            case 4:  res = 0.5*pow(zp,2.0) + 0.5*pow(zp,3.0) + 0.4583333333333333*pow(zp,4.0) + 0.41666666666666663*pow(zp,5.0) + 0.3805555555555556*pow(zp,6.0);
                break;
            case 5:  res = 0.16666666666666666*pow(pi,2.0) + (-1. + llzp)*zp + (-0.25 + 0.5*llzp)*pow(zp,2.0) + (-0.1111111111111111 + 0.3333333333333333*llzp)*pow(zp,3.0) + (-0.0625 + 0.25*llzp)*pow(zp,4.0) + (-0.04 + 0.2*llzp)*pow(zp,5.0) + (-0.027777777777777776 + 0.16666666666666666*llzp)*pow(zp,6.0);
                break;
            case 6:  res = 0.5*pow(ll2,2.0) - 1.*ll2*llzp - 0.08333333333333333*pow(pi,2.0) + 0.5*zp + 0.0625*pow(zp,2.0) + 0.013888888888888888*pow(zp,3.0) + 0.00390625*pow(zp,4.0) + 0.00125*pow(zp,5.0) + 0.00043402777777777775*pow(zp,6.0);
                break;
            case 7:  res = -0.16666666666666666*pow(pi,2.0) + zp + 0.25*pow(zp,2.0) + 0.1111111111111111*pow(zp,3.0) + 0.0625*pow(zp,4.0) + 0.04*pow(zp,5.0) + 0.027777777777777776*pow(zp,6.0);
                break;
            case 8:  res = 0.5*pow(llzp,2.0);
                break;
                
            default: {cout << "chaplin error: index in HPL2ar1 out of bounds"
                << " : should be 0-8 and it is " << index << endl;
                exit(EXIT_FAILURE);
            }
        }
        //res = chaplin::correct_imaginary_parts(n1,xx,res);
        
        if (imaginary_part_was_modified and imaginary_part_should_be_zero_in_hpl2(index,x)) {
            res = complex<double> (res.real() , 0.0);
        }
        
        
        
        
        
        return res;
    }
    
    
    
    
    
    complex<double> HPL2ar0(int index,const complex<double> & xx){
        
        
        complex<double> ll2 ( log(2.0) ,  0.0);
        complex<double> ll2_sq ( pow( log(2.0) , 2.0 )  ,  0.0);
        complex<double> pi_sq (consts::pi_square ,0.0 );
        
        complex<double> x=xx;
        
        // adding i*epsilon
        bool imaginary_part_was_modified = false;
        if ( x.imag() == 0.0 ) {
            x = x + complex<double> (0.0, 1e-60);
            imaginary_part_was_modified = true;
        }
        complex<double> zp = x;
        complex<double> res(0.0,0.0);
        complex<double> llx = log(x);
        
        switch (index) {
            case 0: res = (pow(zp,2.))/2.0
                - (pow(zp,3.))/2.0
                + (11.0*pow(zp,4.))/24.0
                - (5.0*pow(zp,5.))/12.0
                + (137.0*pow(zp,6.))/360.0
                - (7.0*pow(zp,7.))/20.0
                + (363.0*pow(zp,8.))/1120.0
                - (761.0*pow(zp,9.))/2520.0
                + (7129.0*pow(zp,10.))/25200.0
                ;
                break;
            case 1: res = -x
                + (pow(zp,2.))/4.0
                - (pow(zp,3.))/9.0
                + (pow(zp,4.))/16.0
                - (pow(zp,5.))/25.0
                + (pow(zp,6.))/36.0
                - (pow(zp,7.))/49.0
                + (pow(zp,8.))/64.0
                - (pow(zp,9.))/81.0
                + (pow(zp,10.))/100.0
                + x*llx
                - (pow(zp,2.)*llx)/2.0
                + (pow(zp,3.)*llx)/3.0
                - (pow(zp,4.)*llx)/4.0
                + (pow(zp,5.)*llx)/5.0
                - (pow(zp,6.)*llx)/6.0
                + (pow(zp,7.)*llx)/7.0
                - (pow(zp,8.)*llx)/8.0
                +(pow(zp,9.)*llx)/9.0
                - (pow(zp,10.)*llx)/10.0;
                break;
            case 2: res = (pow(zp,2.))/2.0
                - (pow(zp,3.))/6.0
                + (5.0*pow(zp,4.))/24.0
                - (7.0*pow(zp,5.))/60.0
                + (47.0*pow(zp,6.))/360.0
                - (37.0*pow(zp,7.))/420.0
                + (319.0*pow(zp,8.))/3360.0
                - (533.0*pow(zp,9.))/7560.0
                + (1879.0*pow(zp,10.))/25200.0;
                break;
            case 3: res =  x
                - (pow(zp,2.))/4.0
                + (pow(zp,3.))/9.0
                - (pow(zp,4.))/16.0
                + (pow(zp,5.))/25.0
                - (pow(zp,6.))/36.0
                + (pow(zp,7.))/49.0
                - (pow(zp,8.))/64.0
                + (pow(zp,9.))/81.0
                - (pow(zp,10.))/100.0;
                break;
            case 4: res =  (llx*llx)/2.0;
                break;
            case 5: res =  x
                + (pow(zp,2.))/4.0
                + (pow(zp,3.))/9.0
                + (pow(zp,4.))/16.0
                + (pow(zp,5.))/25.0
                + (pow(zp,6.))/36.0
                + (pow(zp,7.))/49.0
                + (pow(zp,8.))/64.0
                + (pow(zp,9.))/81.0
                + (pow(zp,10.))/100.0;
                break;
            case 6: res =  (pow(zp,2.))/2.0
                + (pow(zp,3.))/6.0
                + (5.0*pow(zp,4.))/24.0
                + (7.0*pow(zp,5.))/60.0
                + (47.0*pow(zp,6.))/360.0
                + (37.0*pow(zp,7.))/420.0
                + (319.0*pow(zp,8.))/3360.0
                + (533.0*pow(zp,9.))/7560.0
                + (1879.0*pow(zp,10.))/25200.0;
                break;
            case 7: res =  -x
                - (pow(zp,2.))/4.0
                - (pow(zp,3.))/9.0
                - (pow(zp,4.))/16.0
                - (pow(zp,5.))/25.0
                - (pow(zp,6.))/36.0
                - (pow(zp,7.))/49.0
                - (pow(zp,8.))/64.0
                - (pow(zp,9.))/81.0
                - (pow(zp,10.))/100.0
                + x*llx
                + (pow(zp,2.)*llx)/2.0
                + (pow(zp,3.)*llx)/3.0
                + (pow(zp,4.)*llx)/4.0
                + (pow(zp,5.)*llx)/5.0
                + (pow(zp,6.)*llx)/6.0
                + (pow(zp,7.)*llx)/7.0
                + (pow(zp,8.)*llx)/8.0
                +(pow(zp,9.)*llx)/9.0
                + (pow(zp,10.)*llx)/10.0;
                break;
            case 8: res = (pow(zp,2.))/2.0
                + (pow(zp,3.))/2.0
                + (11.0*pow(zp,4.))/24.0
                + (5.0*pow(zp,5.))/12.0
                + (137.0*pow(zp,6.))/360.0
                + (7.0*pow(zp,7.))/20.0
                + (363.0*pow(zp,8.))/1120.0
                + (761.0*pow(zp,9.))/2520.0
                + (7129.0*pow(zp,10.))/25200.0 ;
                break;
            default: {cout << "chaplin error: index in HPL2ar1 out of bounds"
                << " : should be 0-8 and it is " << index << endl;
                exit(EXIT_FAILURE);
            }
        }
        //res = chaplin::correct_imaginary_parts(n1,xx,res);
        
        if (imaginary_part_was_modified and imaginary_part_should_be_zero_in_hpl2(index,x)) {
            res = complex<double> (res.real() , 0.0);
        }
        
        
        
        
        return res;
    }
    
    complex<double> HPL2arm1(int index,const complex<double> & xx){
        
        complex<double> ll2 ( log(2.0) ,  0.0);
        complex<double> ll2_sq ( pow( log(2.0) , 2.0 )  ,  0.0);
        complex<double> pi_sq (consts::pi_square ,0.0 );
        const complex<double> pi(consts::Pi,0.0);
        const complex<double> myi(0.0,1.0);
        complex<double> x=xx;
        
        // adding i*epsilon
        bool imaginary_part_was_modified = false;
        if ( x.imag() == 0.0 ) {
            x = x + complex<double> (0.0, 1e-60);
            imaginary_part_was_modified = true;
        }
        
        complex<double> zp = 1.0 + x;
        complex<double> llzp = log(zp);
        complex<double> res(0.0,0.0);
        complex<double> szp = mysign(zp);
        
        switch (index) {
            case 0:  res = 0.5*pow(llzp,2.0);
                break;
            case 1:  res = 0.16666666666666666*pow(pi,2.0) + llzp*myi*pi*szp - zp - 0.25*pow(zp,2.0) - 0.1111111111111111*pow(zp,3.0) - 0.0625*pow(zp,4.0) - 0.04*pow(zp,5.0) - 0.027777777777777776*pow(zp,6.0) - 0.02040816326530612*pow(zp,7.0) - 0.015625*pow(zp,8.0) - 0.012345679012345678*pow(zp,9.0);
                break;
            case 2:  res = 0.5*pow(ll2,2.0) - ll2*llzp - 0.08333333333333333*pow(pi,2.0) + 0.5*zp + 0.0625*pow(zp,2.0) + 0.013888888888888888*pow(zp,3.0) + 0.00390625*pow(zp,4.0) + 0.00125*pow(zp,5.0) + 0.00043402777777777775*pow(zp,6.0) + 0.00015943877551020407*pow(zp,7.0) + 0.00006103515625*pow(zp,8.0) + 0.000024112654320987653*pow(zp,9.0);
                break;
            case 3:  res = -0.16666666666666666*pow(pi,2.0) + (1. - llzp)*zp + (0.25 - 0.5*llzp)*pow(zp,2.0) + (0.1111111111111111 - 0.3333333333333333*llzp)*pow(zp,3.0) + (0.0625 - 0.25*llzp)*pow(zp,4.0) + (0.04 - 0.2*llzp)*pow(zp,5.0) + (0.027777777777777776 - 0.16666666666666666*llzp)*pow(zp,6.0) + (0.02040816326530612 - 0.14285714285714285*llzp)*pow(zp,7.0) + (0.015625 - 0.125*llzp)*pow(zp,8.0) + (0.012345679012345678 - 0.1111111111111111*llzp)*pow(zp,9.0);
                break;
            case 4:  res = -1.0*0.5*pow(pi,2.0) - myi*pi*szp*zp + (0.5 - 0.5*myi*pi*szp)*pow(zp,2.0) + (0.5 - 0.3333333333333333*myi*pi*szp)*pow(zp,3.0) + (0.4583333333333333 - 0.25*myi*pi*szp)*pow(zp,4.0) + (0.41666666666666663 - 0.2*myi*pi*szp)*pow(zp,5.0) + (0.3805555555555556 - 0.16666666666666666*myi*pi*szp)*pow(zp,6.0) + (0.35000000000000003 - 0.14285714285714285*myi*pi*szp)*pow(zp,7.0) + (0.32410714285714287 - 0.125*myi*pi*szp)*pow(zp,8.0) + (0.30198412698412697 - 0.1111111111111111*myi*pi*szp)*pow(zp,9.0);
                break;
            case 5:  res = -1.0*0.08333333333333333*pow(pi,2.0) + ll2*zp + (-0.25 + 0.5*ll2)*pow(zp,2.0) + (-0.20833333333333331 + 0.3333333333333333*ll2)*pow(zp,3.0) + (-0.16666666666666666 + 0.25*ll2)*pow(zp,4.0) + (-0.13645833333333332 + 0.2*ll2)*pow(zp,5.0) + (-0.11475694444444445 + 0.16666666666666666*ll2)*pow(zp,6.0) + (-0.09873511904761906 + 0.14285714285714285*ll2)*pow(zp,7.0) + (-0.08653273809523811 + 0.125*ll2)*pow(zp,8.0) + (-0.07697224289021165 + 0.1111111111111111*ll2)*pow(zp,9.0);
                break;
            case 6:  res = -1.0*0.5*pow(ll2,2.0) + 0.08333333333333333*pow(pi,2.0) + (-0.5 + 0.5*llzp)*zp + (-0.0625 + 0.125*llzp)*pow(zp,2.0) + (-0.013888888888888888 + 0.041666666666666664*llzp)*pow(zp,3.0) + (-0.00390625 + 0.015625*llzp)*pow(zp,4.0) + (-0.00125 + 0.00625*llzp)*pow(zp,5.0) + (-0.00043402777777777775 + 0.0026041666666666665*llzp)*pow(zp,6.0) + (-0.00015943877551020407 + 0.0011160714285714285*llzp)*pow(zp,7.0) + (-0.00006103515625 + 0.00048828125*llzp)*pow(zp,8.0) + (-0.000024112654320987653 + 0.00021701388888888888*llzp)*pow(zp,9.0);
                break;
            case 7:  res = 0.08333333333333333*pow(pi,2.0) - ll2*myi*pi*szp + 0.5*myi*pi*szp*zp + (-0.25 + 0.125*myi*pi*szp)*pow(zp,2.0) + (-0.16666666666666666 + 0.041666666666666664*myi*pi*szp)*pow(zp,3.0) + (-0.10416666666666666 + 0.015625*myi*pi*szp)*pow(zp,4.0) + (-0.06666666666666667 + 0.00625*myi*pi*szp)*pow(zp,5.0) + (-0.044444444444444446 + 0.0026041666666666665*myi*pi*szp)*pow(zp,6.0) + (-0.030952380952380953 + 0.0011160714285714285*myi*pi*szp)*pow(zp,7.0) + (-0.0224702380952381 + 0.00048828125*myi*pi*szp)*pow(zp,8.0) + (-0.016931216931216932 + 0.00021701388888888888*myi*pi*szp)*pow(zp,9.0);
                break;
            case 8:  res = 0.5*pow(ll2,2.0) - 0.5*ll2*zp + (0.125 - 0.125*ll2)*pow(zp,2.0) + (0.0625 - 0.041666666666666664*ll2)*pow(zp,3.0) + (0.028645833333333332 - 0.015625*ll2)*pow(zp,4.0) + (0.013020833333333332 - 0.00625*ll2)*pow(zp,5.0) + (0.005946180555555556 - 0.0026041666666666665*ll2)*pow(zp,6.0) + (0.0027343750000000003 - 0.0011160714285714285*ll2)*pow(zp,7.0) + (0.0012660435267857143 - 0.00048828125*ll2)*pow(zp,8.0) + (0.000589812748015873 - 0.00021701388888888888*ll2)*pow(zp,9.0);
                break;
            default: {cout << "chaplin error: index in HPL2arm1 out of bounds"
                << " : should be 0-8 and it is " << index << endl;
                exit(EXIT_FAILURE);
            }
        }
        //res = chaplin::correct_imaginary_parts(n1,xx,res);
        
        if (imaginary_part_was_modified and imaginary_part_should_be_zero_in_hpl2(index,x)) {
            res = complex<double> (res.real() , 0.0);
        }
        
        
        
        
        return res;
    }
    
    complex<double> HPL2else(int index,const complex<double> & xx){
        
        complex<double> ll2 ( log(2.0) ,  0.0);
        complex<double> ll2_sq ( pow( log(2.0) , 2.0 )  ,  0.0);
        complex<double> pi_sq (consts::pi_square ,0.0 );
        const complex<double> pi(consts::Pi,0.0);
        const complex<double> myi(0.0,1.0);
        complex<double> x=xx;
        
        // adding i*epsilon
        bool imaginary_part_was_modified = false;
        if ( x.imag() == 0.0 ) {
            x = x + complex<double> (0.0, 1e-60);
            imaginary_part_was_modified = true;
        }
        
        //complex<double> zp = 1.0 + x;
        //complex<double> llzp = log(zp);
        complex<double> res(0.0,0.0);
        complex<double> llx = log(x);
        complex<double> ll1x = log(1.+x);
        complex<double> ll1mx = log(1.-x);
        //complex<double> szp = mysign(zp);
        
        
        
        
        switch (index) {
            case 0:  res = 0.5*pow(ll1x,2.0);
                break;
            case 1:  res = ll1x*llx + basis2x2(x);
                break;
            case 2:  res = -(ll1mx*ll1x) + ll1mx*ll2 - 0.5*pow(ll2,2.0) + 0.08333333333333333*pow(pi,2.0) - basis2x3(x);
                break;
            case 3:  res = -basis2x2(x);
                break;
            case 4:  res = 0.5*pow(llx,2.0);
                break;
            case 5:  res = basis2x1(x);
                break;
            case 6:  res = -(ll1mx*ll2) + 0.5*pow(ll2,2.0) - 0.08333333333333333*pow(pi,2.0) + basis2x3(x);
                break;
            case 7:  res = -(ll1mx*llx) - basis2x1(x);
                break;
            case 8:  res = 0.5*pow(ll1mx,2.0);
                break;
            default: {cout << "chaplin error: index in HPL2else out of bounds"
                << " : should be 0-8 and it is " << index << endl;
                exit(EXIT_FAILURE);
            }
        }
        //res = chaplin::correct_imaginary_parts(n1,xx,res);
        
        if (imaginary_part_was_modified and imaginary_part_should_be_zero_in_hpl2(index,x)) {
            res = complex<double> (res.real() , 0.0);
        }
        
        
        
        
        return res;
    }
    


}