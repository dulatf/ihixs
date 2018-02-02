#include "higgs_eft.h"


using namespace std;

namespace HEFT {

    double one_minus_z(const double& z, const double& L)
    {
        return 1.-z;
    }
    double one(const double& z, const double& L)
    {
        return 1.;
    }

    double intpow(const double& x,int m){
        double res=1.0;
        for (int i=0;i<m;i++){
            res *= x;
        }
        return res;
    }
    
    
    
    
    complex<double> operator*(int i,const complex<double>& c){return double(i)*c;}
    complex<double> operator*(const complex<double>& c,int i){return i*c;}
    
    void check_imaginary_part(const complex<double> c, const char* func_name)
    {
        if (abs(imag(c))>1e-3)
        {
            cout<<"\n nnlo regular term found having imaginary part! "
            <<c
            <<" coming from function: "<<func_name
            <<endl;
        }
    }
    
    
}




