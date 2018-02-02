// Copyright 2018 ihixs team

#ifndef SRC_CORE_INPUT_PARAMETERS_H_
#define SRC_CORE_INPUT_PARAMETERS_H_

#include <string>
#include "src/tools/luminosity.h"
#include "src/core/luminosity_integrals.h"
#include "src/tools/user_interface.h"
#include "src/higgs/effective_theory/wilson_coefficients.h"
#include "src/tools/as_series.h"
#include "src/models/model.h"
#include "src/higgs/electroweak_corrections/gluon_fusion_ew_coefficients.h"

class InputParameters{
 public:
    void Configure(const UserInterface& UI);
    double muf() {return _muf;}
    string InputInformation() {return _input_information;}

 public:
    CModel _model;
    WilsonCoefficient _wc;
    AsSeries _ew_lambda_series;
    double _prefactor;
    double _top_only_LO_coefficient;
    double _tau;
    double _mur;
    double _muf;
    double _as_over_pi;
    double _as_at_mz;
    double _log_muf_over_mt_sq;
    double _log_muf_over_mh_sq;
    double _log_mur_over_muf_sq;
    double _log_muf_over_mt_sq_for_WC;
    int _int_qcd_perturbative_order;
    int _int_qcd_perturbative_order_for_model_evolution;
    bool _with_exact_qcd_corrections;
    bool _with_ew_corrections;
    Luminosity* _lumi;
    VegasArguments _VI;

 private:
    int DetermineQCDPerturbativeOrder(const string& order);
    void SetUpWilsonCoefficient(const UserInterface& UI);
    void ComputeEwLambdaSeries();
    void CreateInputInformation();
    void ComputeTopOnlyLOCoefficient(const UserInterface& UI);

 private:
    string _input_information;
    double _lambda;
};


#endif  // SRC_CORE_INPUT_PARAMETERS_H_
