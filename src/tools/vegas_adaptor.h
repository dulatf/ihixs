

#ifndef VEGAS_ADAPTOR_H
#define VEGAS_ADAPTOR_H

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
using namespace std;
//#include "user_interface.h"
//#include "thehatch.h"
//#include "bin.h"



typedef int (*pointer_to_Integrand)(const int *ndim, const double xx[],
const int *ncomp, double ff[],void * therun, double* weight, int* iteration_number);

//  ---------------------------------------------------------------------------
//
//  classes for easy vegas integration
//
//  usage:  subclass CoolInt to overload the function
//          double evaluateIntegral(const double xx[]);
//          Then use by constructing an object of the derived class
//          setting it up potentially
//          and then call_vegas();
//
//

class CoolInt;

#ifndef POINTER_TO_COOL_INTEGRAL
#define POINTER_TO_COOL_INTEGRAL
//CoolInt * ptr_to_cool_int;
#endif

int cool_integral(const int *ndim, const double xx[],
                  const int *ncomp, double ff[],void * therun,
                  double* weight);


class CoolInt{
public: //data
    

public: //functions
    // default constructor setting default 
    CoolInt();
    ~CoolInt(){};
    
    void setParams(int number_of_dims,double epsrel,double epsabs,
                   int mineval,int maxeval,int nstart,int nincrease);
    void setParams(int number_of_dims,double epsrel,double epsabs);
    
    virtual vector<double> evaluateIntegral(const double xx[]);
    void call_vegas();
    void call_cuhre();
    double result(){return _central_value[0];}
    double error(){return _vegas_error[0];}
    double central_value(){return _central_value[0];}
    double mc_error_of_central_value(){return _vegas_error[0];}

    
    double GiveEpsRel(){return _epsrel*_rel_accuracy_multiplier;}
    
    int Dimension(){return _number_of_dims;}
    //double give_res_component(unsigned i);
    //double give_err_component(unsigned i);

    vector<double> GiveResultVector();
    vector<double> GiveErrorVector();
    
    
    void set_number_of_xs_values_we_keep(unsigned);
    
    double _vegas_weight;
    
    
protected: //functions
    

    void set_dimensions(int dim){_number_of_dims=dim;}
    
    void set_rel_accuracy_multiplier(double x){_rel_accuracy_multiplier=x;}
    

    

protected://data
    int _nbatch,_verbose;
    double _epsrel,_epsabs;
    int _mineval,_maxeval,_nstart,_nincrease;
    
    double _rel_accuracy_multiplier;
    int _number_of_components;
    
private:
    void check_number_of_components();

private: //data
    int _number_of_dims;
    double _central_value[120];
    double _vegas_error[120];
    double _prob[120];
    int _neval,_fail;
    int _gridno;
    int _seed;
    int _spin[1];

};




#endif



