
#include <math.h>       /* round*/
#include <complex>
#include "constants.h" // Pi
#include <iostream>
using namespace std;


//: implementation from K.S. Kölbig,
//: "Programs for computing the logarithm of the gamma function, and the digamma function, for complex argument"
//:
//: but note the wrong inversion formula in the case Re(z)<0. The correct inversion formula is the one implemented here. Checked against the wpsipg.F from CERNlib (agreement to 16 digits). Also checked against Mathematica (that somehow doesn't want to provide more than 6 digits for Re(z) not an integer).

complex<double> psi(const complex<double>& z) {
    double x = z.real();
    if (x<0) return psi(1.-z) - consts::Pi / tan(consts::Pi * z);
    else if (0 <= x and x <7) {
        int n = 7 - floor(x) ;
        //cout << " n = "<< n << endl;
        complex<double> res(0.0,0.0);
        for (int v = 0 ; v < n ; v++) {
            //cout << "adding " << 1./(z + double(v)) << " -->res = "<< res << endl;
            res = res - 1./(z + double(v));
        }
        res = res + psi( z + double(n) );
        return res;
    }
    else if (x >=7) {
        complex<double> res = log(z) - 1./2./z;
        //: Bernoulli numbers from mathematica: BernoulliB[2*k]/(2*k)
        //: with k starting at 1 (so B2k[0] = BernoulliB[2]/2)
        double B2k[15]={0.08333333333333333,-0.008333333333333333,0.003968253968253968,
            -0.004166666666666667,0.007575757575757576,-0.02109279609279609,0.08333333333333333,
            -0.4432598039215686,3.05395433027012,-26.45621212121212,281.4601449275362,
            -3607.510546398046,54827.58333333333,-974936.8238505747,2.005269579668808e7};
        for (int k = 1 ; k < 16 ; k++) {
            res = res - B2k[k-1] * pow(z,-2.*k);
        }
        return res;
    }
   std::cout << "psi(" << z << "): this should not have been reached, invalid input?\n";
   return 0; 
}




