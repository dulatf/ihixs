#include "psi.h"
#include <iostream>
using namespace std;

extern "C" {
        complex<double> wpsipg_(complex<double> *z, int *k); 
}
complex<double> PsiFinder::get(const complex<double>& zz){
    
    complex<double> myz{zz};
    int myk = 0;
    complex<double> cres = wpsipg_(&myz,&myk);
    return cres;

}

