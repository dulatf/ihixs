#include "luminosity_integrals.h"
#include <fstream>
//#include <io.h>   // For access().
#include <sys/types.h>  // For stat().
#include <sys/stat.h>   // For stat().
#include <string>
using namespace std;

void VegasArguments::Configure(const UserInterface& UI){
    epsrel = UI.giveDouble("epsrel");
    epsabs = UI.giveDouble("epsabs");
    mineval = UI.giveInt("mineval");
    maxeval = UI.giveInt("maxeval");
    nstart = UI.giveInt("nstart");
    nincrease = UI.giveInt("nincrease");
    verbose = UI.giveInt("cuba_verbose");}

void LuminosityIntegral::set_initial_flavors(){
    bool found=false;
    if (_channel=="gg") {
        lumi_->clear_pairs();//if you don't clear you will pay the price in factors
        lumi_->addPair(QCD::g, QCD::g);
        found=true;
    }
    if (_channel=="qg") {
        lumi_->clear_pairs();
        for (int i=-5;i<6;i++)
        {
            if (i!=0)
            {
                //:  q g + qbar g
                lumi_->addPair(i, QCD::g);
                //:  g q + g qbar
                lumi_->addPair(QCD::g, i);
            }
        }
        found=true;
    }
    if (_channel=="qqbar"){
        lumi_->clear_pairs();
        for (int i=-5;i<6;i++)
        {
            if (i!=0)
            {
                //:  q qbar + qbar q
                lumi_->addPair(i, -i);
            }
        }
        found=true;
    }
    if (_channel=="qq"){
        lumi_->clear_pairs();
        for (int i=-5;i<6;i++)
        {
            if (i!=0)
            {
                //:  q q
                lumi_->addPair(i, i);
            }
        }
        found=true;
    }
    if (_channel=="q1q2"){
        lumi_->clear_pairs();
        for (int i=-5;i<6;i++)
        {
            for (int j=-5;j<6;j++)
            {
                if (i!=0 and j!=0 and i!=j and i!=-j)
                {
                    //:  q1 q2
                    lumi_->addPair(i, j);
                }
            }
        }
        found=true;
    }
    if (not(found)){
        cout<<"[LuminosityIntegrals]: initial state "<<_channel<<" not defined ??"<<endl;
        exit(1);
    }
    
    
}


void LuminosityIntegral::Configure(Luminosity* lumi,
               const double& tau,
               const VegasArguments& VA,
               const double& muf,
               const double& log_muf_over_mh_sq)
{   lumi_=lumi;
    tau_ = tau;
    set_initial_flavors();
    _epsrel = VA.epsrel;
    _epsabs = VA.epsabs;
    _mineval = VA.mineval;
    _maxeval = VA.maxeval;
    _nstart = VA.nstart;
    _nincrease = VA.nincrease;
    _verbose = VA.verbose;
    _muF = muf;
    _Lf = log_muf_over_mh_sq;
    _number_of_components = lumi_->Size();
    
}

void LuminosityIntegral::PerformIntegration(){
    
        PerformIntegrationByVegas();
    
    
}


void LuminosityIntegral::PerformIntegrationByVegas(){
    if (Dimension()>1 and Dimension()<4 and _integration_method!="vegas"){
        call_cuhre();
    }
    else{
        call_vegas();
    }
    _vegas_result = ResultPair( result(), error()  );
}



//--------------------------------------------------------
vector<double> LuminosityIntegralAbsolutelyGenericDelta::evaluateIntegral(const double* xx)
{
    const double x1= xx[0];
    const double measure = 1./x1;
    const double res= measure*lumi_->give(x1,tau_/x1,_muF);
    
    return vector<double>(1,res);
}



vector<double> LuminosityIntegralAbsolutelyGenericPlus::evaluateIntegral(const double* xx)
{
    const double x1= xx[0];
    const double z = xx[1];
    const double measure = 1./x1;
   
    const double res=  measure
    *( lumi_->give(x1,tau_/x1/z,_muF) - lumi_->give(x1,tau_/x1,_muF) )
    / (1.-z)
    * pow(log(1.-z), _log_power);

    
    
    return vector<double>(1,res);
}


vector<double> LuminosityIntegralAbsolutelyGenericReg::evaluateIntegral(const double* xx)
{
    const double x1= xx[0];
    const double z = xx[1];
    double measure = 1./x1;
    
    compute_matrix_element_Lf_coefficients(xx);
    double sum=0.0;
    for (int i=0;i<_me_lf_coeffs.size();i++) sum += _me_lf_coeffs[i]*pow(_Lf,i);
    
    
     double res=  measure *lumi_->give(x1,tau_/x1/z,_muF)
                        * sum;
    
    
    
    return vector<double>(1,res);
}






