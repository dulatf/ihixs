
#include "vegas_adaptor.h"
//#include "mincuba.h"
#include "cuba.h"
#include "math.h"
#include "sstream"



int Integrand(const int *ndim, const double xx[],
                     const int *ncomp, double ff[],void * therun, double* weight);




vector<double> CoolInt::evaluateIntegral(const double xx[])
{
    cout<<"CoolInt: error: calling default integrand"<<endl;
    exit(0);
    double lambda = xx[0];
    return vector<double>(1,pow(lambda,3.0));
    
}


CoolInt::CoolInt()
{
    _gridno=0;
    _seed=0;//:0: Sobol,  >0 Ranlux
    _nbatch=10;
    _mineval=10000;
    _maxeval=50000000;
    _nstart=5000;
    _nincrease=1000;
    _verbose = 0;
    _number_of_dims = 1;
    _epsrel = 1e-3;
    _epsabs = 0.0;
    _number_of_components=1;
    //cout<<"\n**Setting _ncomp "<<_number_of_components<<endl;

    _rel_accuracy_multiplier=1.0;
    _spin[0]= -1;// see cuba doc for the meaning of this
}


void CoolInt::setParams(int number_of_dims,double epsrel,double epsabs,
                        int mineval,int maxeval,int nstart,int nincrease)

{
    _gridno=0;
    _seed=0;//:0: Sobol,  >0 Ranlux
    _nbatch=10;
    //_verbose = 0;
    _number_of_dims = number_of_dims;
    _number_of_components=1;

    _epsabs = epsabs;
    _epsrel = epsrel;
    _mineval = mineval;
    _maxeval = maxeval;
    _nstart = nstart;
    _nincrease = nincrease;
}

void CoolInt::setParams(int number_of_dims,double epsrel,double epsabs)
{
    _gridno=0;
    _seed=0;//:0: Sobol,  >0 Ranlux
    _nbatch=10;
    _mineval=10000;
    _maxeval=50000000;
    _nstart=5000;
    _nincrease=1000;
    //_verbose = 0;
    _number_of_dims = number_of_dims;
    _number_of_components=1;

    _epsrel = epsrel;
    _epsabs = epsabs;
}


void CoolInt::call_vegas()
{
    check_number_of_components();
    //:cuba 4.2 interface
    Vegas(_number_of_dims,
          _number_of_components,
          integrand_t(&cool_integral),
          this,
          1,// nvec (see cuba documentation for what this does)
          GiveEpsRel(),
          _epsabs,
          _verbose,
          _seed,
          _mineval,
          _maxeval,
          _nstart,
          _nincrease,
          _nbatch,
          _gridno,
          NULL,//grid_file_name.c_str(),
          _spin,//,void *spin (see cuba documentation for what this does)
          &_neval,
          &_fail,
          _central_value,
          _vegas_error,
          _prob);
    
    
}

void CoolInt::call_cuhre(){
    check_number_of_components();

     int cuhre_key=9;
    // double cuhre_central[1];
    // double cuhre_error[1];
    // double cuhre_prob[1];
     int _nregions;
    double cuhre_mineval=1000;
    
   

     //:cuba 4.2 interface
    Cuhre(_number_of_dims,
               _number_of_components,
               integrand_t(&cool_integral),
               this,
               1,// nvec (see cuba documentation for what this does)
               GiveEpsRel(),
               _epsabs,
               _verbose,
               cuhre_mineval,
               _maxeval,
               cuhre_key,
               NULL,//grid_file_name.c_str(),
               _spin,//,void *spin (see cuba documentation for what this does)
               &_nregions,
               &_neval,
               &_fail,
               _central_value,
               _vegas_error,
               _prob);
     
     
  
}

void CoolInt::check_number_of_components(){
    if (_number_of_components>120){
        cout<<"Error in vegas interface (CoolInt): number of components "<<_number_of_components
            <<" too large. Maximum allowed is 120. Change the magic number in vegas_adaptor.h for more."
            <<endl;
        exit(0);
    }
}



vector<double> CoolInt::GiveResultVector(){
    vector<double> res;
    for (int i=0;i<_number_of_components;i++){
        res.push_back(_central_value[i]);
    }
    return res;
}

vector<double> CoolInt::GiveErrorVector(){
    vector<double> res;
    for (int i=0;i<_number_of_components;i++){
        res.push_back(_vegas_error[i]);
    }
    return res;
}


int cool_integral(const int *ndim, const double xx[],
                  const int *ncomp, double ff[],void * theclass,
                  double* weight_from_vegas)
{
    CoolInt* ptr_to_class = static_cast<CoolInt*>(theclass);
    ptr_to_class->_vegas_weight = *weight_from_vegas;
    
    vector<double> integrand_values = ptr_to_class->evaluateIntegral(xx);
    
    for (int i=0;i<integrand_values.size();i++){
        //cout<<"computing component no. "<<i<<endl;
        ff[i] = integrand_values[i];
        }
    if (integrand_values.size()<*ncomp){
        for (int i=integrand_values.size();i<*ncomp;i++)
            ff[i]=0.0;
        
        cout<<"\nwarning: ncomp > integrand_values.size()"
            <<"we set the components of the integrand from "<<integrand_values.size()
            <<" to "<<*ncomp<<" to zero "<<endl;
    }
    
    if (*ncomp<integrand_values.size()){
        cout<<"Error in CoolInt: integran_values = "<<integrand_values.size()<<" > ncomp = "<<*ncomp<<endl;
        exit(0);
    }
    
    return(0);
}










