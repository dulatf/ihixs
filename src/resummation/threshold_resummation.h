#ifndef THRESHOLD_RESUMMATION_H
#define THRESHOLD_RESUMMATION_H

#include "constants.h"
#include <vector>
#include "as_series.h"
#include "user_interface.h"
#include "input_parameters.h"
#include <complex>

#include "gauss_grid.h"
using namespace std;


class ThresholdResummation{
public:
  ThresholdResummation(const UserInterface& UI, const string &gridDirectory, InputParametersForThresRes *pInput);
  ~ThresholdResummation();
  double ResummationCorrection(unsigned int logOrder, unsigned int matchingOrder, bool pi2Resummation);
private:
  void Initialize(const UserInterface& UI, const string &gridDirectory);
  string generate_grid_name(const UserInterface &UI);


  double Ch(unsigned int ord, bool piSquareResummed);
  complex<double> MatchingCoefficient(const complex<double> &n, unsigned int pertOrd, unsigned int logOrd);
  complex<double> SudakovExponential(const complex<double> &n, unsigned int lOrd, bool pi2Resummation);

  GaussGrid *_theGrid;
  InputParametersForThresRes *_input;
    
    string _resummation_type;
};



#endif
