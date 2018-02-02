

#ifndef LUMINOSITY_INTEGRALS_H
#define LUMINOSITY_INTEGRALS_H

#include "vegas_adaptor.h"
#include "luminosity.h"
#include "user_interface.h"
#include "constants.h" // for QCD namespace
#include "result_pair.h"

class VegasArguments{
public:
    void Configure(const UserInterface& UI);
    
    double epsrel;
    double epsabs;
    int mineval;
    int maxeval;
    int nstart;
    int nincrease;
    int verbose;
};


class LuminosityIntegral: public CoolInt
{
public:
    virtual vector<double> evaluateIntegral(const double* xx)=0;
    virtual void set_initial_flavors();
    void Configure(Luminosity* lumi,
                   const double& tau,
                   const VegasArguments& VA,
                   const double& muf,
                   const double& log_muf_over_mh_sq);
    double tau(){return tau_;}
    
    
    void PerformIntegration();
    
    ResultPair GiveResult(){return _vegas_result;}
    string Name(){return _name;}
    void set_integration_method(const string& method){_integration_method=method;}
protected:
    Luminosity* lumi_;
    double tau_;
    double _muF;
    double _Lf;
    string _channel;
    
    vector<double> _me_lf_coeffs;

    ResultPair _vegas_result;
    
    string _name;
    
    string _integration_method;
private:
    void PerformIntegrationByVegas();
    
    
    
    
};

//----------------------------------------------------------------


class LuminosityIntegralAbsolutelyGenericDelta: public LuminosityIntegral
{
public:
    vector<double>  evaluateIntegral(const double* xx);
    
};




class LuminosityIntegralAbsolutelyGenericPlus: public LuminosityIntegral
{
public:
    vector<double>  evaluateIntegral(const double* xx);
protected:
    int _log_power;
};

class LuminosityIntegralAbsolutelyGenericReg: public LuminosityIntegral
{
public:
    vector<double>  evaluateIntegral(const double* xx);
    virtual void compute_matrix_element_Lf_coefficients(const double* xx)=0;
    
};








#endif
