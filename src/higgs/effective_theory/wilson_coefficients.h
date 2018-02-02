
#ifndef WILSON_COEFFICIENTS_H
#define WILSON_COEFFICIENTS_H

#include<vector>
#include<string>
#include "as_series.h"
using namespace std;

class WilsonCoefficient{
public:
    void Configure(const double& log_muf_over_mt_sq, const string& scheme);
    AsSeries c() const {return _c;}
    
private:
    vector<double> _w;
    AsSeries _c;
    double _c0,_c1,_c2,_c3;
};

#endif

