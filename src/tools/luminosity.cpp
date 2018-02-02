//#include <string>
//#include <iostream>
//#include <vector>
//#include <math.h>
//#include <stdlib.h>



#include "luminosity.h"
//#include<algorithm>

const double Luminosity::_almost_zero = 1e-12;
/// warning: numerical value of _almost_zero is arbitrary here

Luminosity::Luminosity(const string& gridname)
{
    LHAPDF::setVerbosity(0);

    // gridname example: MSTW2008nnlo68cl
    _pdf = LHAPDF::mkPDFs(gridname);
    // \warning _pdf stores now pointers to all pdf members of this set (e.g. 101 for PDF4LHC15_nnlo_100)
    /// \warning  what happens if the gridname is not matching the pdf names in LHAPDF?
}

Luminosity::Luminosity(const string& gridname,int imember)
{
    LHAPDF::setVerbosity(0);


    
    // gridname example: MSTW2008nnlo68cl
    _pdf.clear();
    _pdf.push_back(LHAPDF::mkPDF( gridname, imember));
    /// \warning  what happens if the gridname is not matching the pdf names in LHAPDF?
    
    //cout<<"succeeded"<<endl;
}

double Luminosity::give(const double& x1,const double& x2,const double& muf)
{
    
    return give_specific(x1,x2,muf,0);
}

vector<double>Luminosity::luminosity_vector(const double& x1,const double& x2,const double& muf){
    vector<double> res;
    for (int i=0;i<_pdf.size();i++){
        res.push_back(give_specific(x1,x2,muf,i));
    }
    return res;
}


double Luminosity::give_specific(const double& x1,const double& x2,const double& muf,int imember)
{
    if (x1>1.0-_almost_zero or x2>1.0-_almost_zero or x1<_almost_zero or x2<_almost_zero)
    {
        return 0.0;
    }
    double res=0.0;
    for (int i=0;i< _pairs.size();i++)
    {
        res +=    _pdf[imember]->xfxQ(_pairs[i].first,x1,muf)
        * _pdf[imember]->xfxQ(_pairs[i].second,x2,muf)
        * _coeff[i];
    }
    return res;
}

double Luminosity::give_single(const double &x, const double &muf, int flavor){
  if(x > 1.0-_almost_zero or x < _almost_zero)
    return 0.0;
  return _pdf[0]->xfxQ(flavor,x,muf);
}

void Luminosity::addPair(int left,int right)
{
    _pairs.push_back(pair<int,int>(left,right));
    _coeff.push_back(1.);
}

void Luminosity::addPair(int left,int right,const double& c)
{
    _pairs.push_back(pair<int,int>(left,right));
    _coeff.push_back(c);
}


#include "LHAPDF/LHAPDF.h"
#include "boost/assign.hpp"
#include <iostream>
#include <fstream>
#include <limits>
using namespace LHAPDF;
using namespace std;
using namespace boost::assign;

double Luminosity::as_at(const double mur){
    const double inf = numeric_limits<double>::infinity();

    LHAPDF::AlphaS_ODE as_ode;
    // As above: order = 0 returns
    // constant value set by
    // as_ode.setAlphaSMZ(double value);
    as_ode.setOrderQCD(4);
    as_ode.setMZ(91);
    as_ode.setAlphaSMZ(alpha_s_at_mz());
    //  as_ode.setMassReference(4.1);
    //  as_ode.setAlphaSReference(0.21);
    as_ode.setQuarkMass(1, 0.0017);
    as_ode.setQuarkMass(2, 0.0041);
    as_ode.setQuarkMass(3, 0.1);
    as_ode.setQuarkMass(4, 1.29);
    as_ode.setQuarkMass(5, 4.1);
    as_ode.setQuarkMass(6, 172.5);
    
    
    const double as_ode_q = as_ode.alphasQ(mur);
    cout << "ODE solution:                  " << setprecision(6) << setw(6)  << ( (as_ode_q > 2) ? inf : as_ode_q )
    << "    num flavs = " << as_ode.numFlavorsQ(mur) << endl;
    cout << mur << " " << as_ode_q << endl;
    return as_ode_q;
}
