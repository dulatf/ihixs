// Copyright 2015 ihixs team
#include "src/core/sigma_term.h"
#include "iomanip"



void SigmaTerm::Evaluate(const InputParameters& input) {
    // initializing luminosity
     _lumi_int->Configure(input._lumi, input._tau, input._VI,
                          input._muf,
                          input._log_muf_over_mh_sq);
    // performing the luminoity integral
    CallVegas();
}



void SigmaTerm::CallVegas() {
    if (not(_evaluated)) {
        //cout << "Computing " << type() << endl;

        _lumi_int->PerformIntegration();

        _result = _result * _lumi_int->GiveResult();
        // the following only concerns the fast run_mode
        _result_vector = _lumi_int->GiveResultVector();
        _error_vector = _lumi_int->GiveErrorVector();
        

        // until here
        _evaluated = true;
    }
    else {
        //cout << "term already evaluated" << endl;
    }
}

ostream& operator<<(ostream& stream, const SigmaTerm& st) {
    stream << setw(28) << left << st._type;
    stream
    << "\t" << st._qcd_result.term_of_order(2)
    << "\t" << st._qcd_result.term_of_order(3)
    << "\t" << st._qcd_result.term_of_order(4)
    << "\t" << st._qcd_result.term_of_order(5)
    << endl;
    return stream;
}


void SigmaTerm::Truncate(int n) {
    _result.Truncate(n);
    _qcd_result.Truncate(n);
    _ew_result.Truncate(n);
}


bool SigmaTerm::IsZero(int porder) {
    if (fabs(_qcd_result.term_of_order(porder).val()) < 1e-14) return true;
    else
        return false;
}

bool SigmaTerm::IsZero() {
    return _result.IsZero();
}















