#include "cppchaplin.h"
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <iostream> /* cout */
#include "src/tools/constants.h"


namespace chaplin{

    complex<double> HPL1at0(int n1);
    complex<double> HPL1at1(int n1);
    complex<double> HPL1atm1(int n1);
    complex<double> HPL1ar1(int n1,const complex<double> & xx);
    complex<double> HPL1ar0(int n1,const complex<double> & xx);
    complex<double> HPL1arm1(int n1,const complex<double> & xx);
    complex<double> HPL1else(int n1,const complex<double> & xx);
    complex<double> correct_imaginary_parts(int n1,
                                            const complex<double>& xx,
                                            const complex<double>& res);
    
    complex<double> HPL(int n1, const complex<double>& x){
        chaplin::check_index(n1);
        
        int region = chaplin::find_region(x);
        
        if ( region == 0 ) {
            return chaplin::HPL1at0(n1);
        }
        else if ( region == 1 ) {
            return chaplin::HPL1at1(n1);
        }
        else if ( region == 2 ) {
            return chaplin::HPL1atm1(n1);
        }
        else if ( region == 3 ) {
            return chaplin::HPL1ar1(n1,x);
        }
        else if ( region == 4 ) {
            return chaplin::HPL1arm1(n1,x);
        }
        else if ( region == 5 ) {
            return chaplin::HPL1ar0(n1,x);
        }
        else return chaplin::HPL1else(n1,x);
    }


    complex<double> correct_imaginary_parts(int n1,
                                            const complex<double>& xx,
                                            const complex<double>& res) {
        
        complex<double> flatres(res.real(),0.0);
        // setting the imaginary part back to zero if it was modified
        if ( xx.imag() == 0 ){
            if ( n1 == 0 and xx.real() > 0) return flatres;
            if ( n1 == 1 and xx.real() < 1) return flatres;
            if ( n1 == -1 and xx.real() > -1) return flatres;
        }
        return res;
    }
    
    


    complex<double> HPL1at0(int n1){
        if (n1==0){
            cout << "chaplin error: call for HPL1(" << n1 << ",0.0)"
            << " which is divergent! "
            << endl;
            exit(EXIT_FAILURE);
        }
        return complex<double>(0.0,0.0);
    }



    complex<double> HPL1at1(int n1){
        if (n1==1){
            cout << "chaplin error: call for HPL1(" << n1 << ",1.0)"
            << " which is divergent! "
            << endl;
            
            exit(EXIT_FAILURE);
        }
        else if (n1 == -1) return complex<double> (log(2.0),0.0);
        else if (n1 ==  0) return complex<double> (0.0,0.0);
        else {
            cout << "chaplin error: call for HPL1(n1,x) with n1 different from {-1,0,1}"
            << "\t n1 = " << n1 << endl;
            exit(EXIT_FAILURE);
        }
    }



    complex<double> HPL1atm1(int n1){
        if (n1 == -1){
            cout << "chaplin error: call for HPL1(" << n1 << ",-1.0)"
            << " which is divergent! "
            << endl;
            
            exit(EXIT_FAILURE);
        }
        else if (n1 == 0) return complex<double> (0.0,consts::Pi);
        else if (n1 == 1) return complex<double> (log(2.0),0.0);
        else {
            cout << "chaplin error: call for HPL1(n1,x) with n1 different from {-1,0,1}"
            << "\t n1 = " << n1 << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    
    
    complex<double> HPL1ar1(int n1, const complex<double>& xx) {
        complex<double> x=xx;
        // adding i*epsilon
        if ( x.imag() == 0.0 ) x = x + complex<double> (0.0, 1e-60);
        complex<double> zp = 1.0-x;
        complex<double> res(0.0,0.0);
        if ( n1 == -1 ) {
            
            
            res =  -((zp)/2.0)
            - pow(zp,2)/8.0
            - pow(zp,3)/24.0
            - pow(zp,4)/64.0
            - pow(zp,5)/160.0
            - pow(zp,6)/384.0
            + log(2.0);
        }
        if ( n1 == 0) {
            
            res =   -zp
            - pow(zp,2)/2.0
            - pow(zp,3)/3.0
            - pow(zp,4)/4.0
            - pow(zp,5)/5.0
            - pow(zp,6)/6.0;
        }
        if ( n1 == 1 ){
            res =  -log(zp);
        }
        
        res = chaplin::correct_imaginary_parts(n1,xx,res);
        return res;
    }
    
    
    
    
    
    
    complex<double> HPL1arm1(int n1, const complex<double>& xx) {
        complex<double> x=xx;
        // adding i*epsilon
        if ( x.imag() == 0.0 ) x = x + complex<double> (0.0, 1e-60);
        complex<double> zp = 1.0+x;
        complex<double> res(0.0,0.0);
        if ( n1 == -1 ) {
            
            
            res =  log( zp );
        }
        if ( n1 == 0) {
            
            
            res = complex<double>(0.0,1.0) * consts::Pi * chaplin::mysign(zp)
            - zp
            - pow(zp,2)/2.0
            - pow(zp,3)/3.0
            - pow(zp,4)/4.0
            - pow(zp,5)/5.0
            - pow(zp,6)/6.0
            - pow(zp,7)/7.0
            - pow(zp,8)/8.0
            - pow(zp,9)/9.0;
        }
        if ( n1 == 1 ){
            //
            res = (zp)/2.0
            + pow(zp,2)/8.0
            + pow(zp,3)/24.0
            + pow(zp,4)/64.0
            + pow(zp,5)/160.0
            + pow(zp,6)/384.0
            + pow(zp,7)/896.0
            + pow(zp,8)/2048.0
            + pow(zp,9)/4608.0
            - log(2.0);
        }
        
        res = chaplin::correct_imaginary_parts(n1,xx,res);
        return res;
    }
    
    
    
    
    complex<double> HPL1ar0(int n1, const complex<double>& xx) {
        complex<double> x=xx;
        // adding i*epsilon
        if ( x.imag() == 0.0 ) x = x + complex<double> (0.0, 1e-60);
        complex<double> res(0.0,0.0);
        if ( n1 == -1 ) {
            
            
            res =  x
            - pow(x,2)/2.0
            + pow(x,3)/3.0
            - pow(x,4)/4.0
            + pow(x,5)/5.0
            - pow(x,6)/6.0
            + pow(x,7)/7.0
            - pow(x,8)/8.0
            + pow(x,9)/9.0
            - pow(x,10)/10.0;
        }
        if ( n1 == 0) {
            res = log(x);
        }
        if ( n1 == 1 ){
            //
            res = x
            + pow(x,2)/2.0
            + pow(x,3)/3.0
            + pow(x,4)/4.0
            + pow(x,5)/5.0
            + pow(x,6)/6.0
            + pow(x,7)/7.0
            + pow(x,8)/8.0
            + pow(x,9)/9.0
            + pow(x,10)/10.0;
        }
        
        res = chaplin::correct_imaginary_parts(n1,xx,res);
        return res;
    }
    
    
    
    complex<double> HPL1else(int n1, const complex<double>& xx) {
        complex<double> x=xx;
        // adding i*epsilon
        if ( x.imag() == 0.0 ) x = x + complex<double> (0.0, 1e-60);
        complex<double> res(0.0,0.0);
        if ( n1 == -1 ) {
            
            
            res =  log( 1.0 + x );
        }
        if ( n1 == 0) {
            res = log(x);
        }
        if ( n1 == 1 ){
            //
            res = - log( 1. - x );
        }
        
        res = chaplin::correct_imaginary_parts(n1,xx,res);
        return res;
    }
    
    

}