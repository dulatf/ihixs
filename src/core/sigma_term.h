// Copyright 2015 ihixs team
#ifndef SRC_CORE_SIGMA_TERM_H_
#define SRC_CORE_SIGMA_TERM_H_

#include "string"
#include "vector"
// using namespace std;
#include "src/tools/luminosity.h"
#include "src/core/luminosity_integrals.h"
#include "src/higgs/effective_theory/wilson_coefficients.h"
// the dependence above is structuraly strange
// should be removed upon refactoring
// or else the SigmaTerm is specific to Higgs production
#include "src/tools/user_interface.h"
#include "src/tools/as_series.h"
#include "src/core/input_parameters.h"

class SigmaTerm{
 public:
    SigmaTerm(
              const string& thetype,
              const AsSeries& val,
              LuminosityIntegral* lumi_int)
    :_type(thetype),
    _result(val),
    _lumi_int(lumi_int),
    _evaluated(false)
        {};
    void CallVegas();
    void Truncate(int order_of_truncation);
    bool IsZero(int porder);
    void Evaluate(const InputParameters& input);
    void ComputeResults(const InputParameters& input);
    void StoreQCDResult(const AsSeries& res) {_qcd_result = res;}
    void StoreEWResult(const AsSeries& res) {_ew_result = res;}
    void StoreEta(const AsSeries& eta) {_eta = eta;}
    void StoreWC(const AsSeries& wc) {_wc = wc;}
    string Name() {return _lumi_int->Name();}
    AsSeries PostVegasResult() {return _result;}
    AsSeries QCDResult() {return _qcd_result;}
    AsSeries EWResult() {return _ew_result;}
    AsSeries Eta() {return _eta;}
    AsSeries WC() {return _wc;}
    LuminosityIntegral* LumiIntegralPtr() {return _lumi_int;}
    friend ostream& operator<<(ostream&, const SigmaTerm&);
    string type() {return _type;}
    bool Evaluated() {return _evaluated;}
    bool IsZero();
    // the following only concerns the fast run_mode
    vector<double> ResultVector() {return _result_vector;}
    vector<double> ErrorVector() {return _error_vector;}
    // until here
 protected:
    string _type;  // delta, plus, reg
    AsSeries _result;
    AsSeries _qcd_result;
    AsSeries _ew_result;
    AsSeries _eta;
    AsSeries _wc;
    LuminosityIntegral* _lumi_int;
    bool _evaluated;
    // the following only concerns the fast run_mode
    vector<double> _result_vector;
    vector<double> _error_vector;
    // until here
};

#endif  // SRC_CORE_SIGMA_TERM_H_"


