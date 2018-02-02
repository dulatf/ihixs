
#ifndef GLUON_FUSION_EW_COEFFICIENTS
#define GLUON_FUSION_EW_COEFFICIENTS

#include "model.h"
#include<complex>
#include<vector>
using namespace std;

struct EwData{
    EwData(const double&m,const double& w){mass=m;deltaew=w;}
    double mass;
    double deltaew;
};

class GluonFusionEWCoefficients
{
public:
    GluonFusionEWCoefficients(const CModel&);
    complex<double> LO(){return NLO_ew_coeff_;}

private:
    complex<double> NLO_ew_coeff_;
    vector<EwData*> ew_data;
private://methods
    vector<double> givecoeff(double x[3],double y[3]);
    
};

#endif

