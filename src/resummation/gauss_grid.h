#ifndef GAUSS_GRID_H
#define GAUSS_GRID_H

#include "constants.h"
#include <vector>
#include "as_series.h"
#include "user_interface.h"
#include <complex>
#include <utility>
#include "luminosity.h"
#include <functional>
using namespace std;


class GaussGrid{
public:
  GaussGrid(const UserInterface& UI,const string &gridPath);
  ~GaussGrid();
  double DoIntegral(const function<complex<double> (complex<double>)>&);
private:
  bool check_grid(const string &gridPath);

  void load_grid(void);
  void store_grid(void);
  void generate_grid(const UserInterface& UI);
  void export_grid(const string& path);
  void generate_subdivisions(void);
  void compute_moments(void);
  complex<double> mellin_integral(const complex<double> &n);
  string _gridPath;

  unsigned int _numContourPoints;
  unsigned int _numIntegralSamples;
  vector<complex<double> > _moments;
  vector<pair<double, double> >_subdivisions;
  Luminosity *_lumi;
  double _scale;
  double _tau;
    bool _verbose;
};

#endif
