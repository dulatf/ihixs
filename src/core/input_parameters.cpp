// Copyright 2018 ihixs team

#include "src/core/input_parameters.h"
#include "src/higgs/exact_qcd_corrections/nlo_exact_matrix_elements.h"


void InputParameters::Configure(const UserInterface& UI) {
    
    _VI.Configure(UI);
    
    _tau = pow(UI.giveDouble("m_higgs"), 2.)/pow(UI.giveInt("Etot"), 2.);
    // chosing whether we do single pdf or all pdfs at the same time
    if (UI.giveBool("with_pdf_error"))
            _lumi = new Luminosity(UI.giveString("pdf_set"));
    else
        _lumi = new Luminosity(UI.giveString("pdf_set"),
                               UI.giveInt("pdf_member"));
    // setting scheme for quarks
    _model.ReadParameters(UI);
    //: 35.0309 = Gf*pi/sqrt(2)/288 with the Gf in pb
    //: Gf = 1.16637*10^{-5} * 0.389379*10^9
    _prefactor = 35.0309;
    // setting up the model and as_pi for the central scale
    _mur = UI.giveDouble("mur");
    _muf = UI.giveDouble("muf");
    // computing log(muf^2/mh^2)
    _log_muf_over_mh_sq = 2. * log(_muf/UI.giveDouble("m_higgs"));
    _with_exact_qcd_corrections = UI.giveBool("with_exact_qcd_corrections");
    _with_ew_corrections = UI.giveBool("with_ew_corrections");
    _int_qcd_perturbative_order = DetermineQCDPerturbativeOrder(
                            UI.giveString("qcd_perturbative_order"));
    _int_qcd_perturbative_order_for_model_evolution =
                            UI.giveInt("qcd_order_evol");
    // note: perturbative order for model evolution: LO is 0 while
    //       _int_qcd_perturbative_order starts at 2
    // new option: fixing the alhpa_s_at_mz externally
    _as_at_mz = _lumi->alpha_s_at_mz();
    if (UI.giveDouble("with_fixed_as_at_mz") != 0.0) {
        _as_at_mz = UI.giveDouble("with_fixed_as_at_mz");
    }
    _model.Configure(
                     _as_at_mz,
                     _mur/UI.giveDouble("m_higgs"),
                     _int_qcd_perturbative_order_for_model_evolution,
                     UI.giveDouble("m_higgs") );
    _as_over_pi = _model.alpha_strong()/consts::Pi;
    _log_muf_over_mt_sq = 2.*log(_muf/_model.top.m());
    _log_mur_over_muf_sq = 2. *  log(_mur/_muf);
    SetUpWilsonCoefficient(UI);
    ComputeTopOnlyLOCoefficient(UI);
    ComputeEwLambdaSeries();
    CreateInputInformation();
}

void InputParameters::SetUpWilsonCoefficient(const UserInterface& UI) {
    // setting up wilson coeffs
    // from NNLO on the wc depends on mu/mt. The top mass for the rest
    // of the computation is run down to mur, but the numerical value of
    // the log should be muf/mt, because that's where the decoupling is
    // performed. This is not related with the a_s running!
    // We need mt at muf
    CModel model_at_muf;
    model_at_muf.ReadParameters(UI);
    // configuring model at muf
    model_at_muf.Configure(
                           _as_at_mz,
                           _muf/UI.giveDouble("m_higgs"),
                           _int_qcd_perturbative_order_for_model_evolution,
                           UI.giveDouble("m_higgs"));
    // setting log(muf^2/mt(mu_f)^2)
    // note: if the OS sheme is used then mt(muf)=mt
    _log_muf_over_mt_sq_for_WC = 2.*log(_muf/model_at_muf.top.m());
    // configure WC with the above log
    _wc.Configure(_log_muf_over_mt_sq_for_WC, _model.top.scheme());
}

void InputParameters::ComputeTopOnlyLOCoefficient(const UserInterface& UI) {
    _model.bottom.set_Y(0.0);
    _model.charm.set_Y(0.0);
    // we compute the LO exact top only matrix element
    h_exact::ExactRegularCalculator top_only(&_model);
    _top_only_LO_coefficient = top_only.born_squared();
    // full exact with bottom and charm
    // setting bottom and charm's Y back to what they
    // were (user defined, defaults to 1.0).
    _model.top.set_Y(UI.giveDouble("y_top"));
    _model.bottom.set_Y(UI.giveDouble("y_bot"));
    _model.charm.set_Y(UI.giveDouble("y_charm"));
}

void InputParameters::ComputeEwLambdaSeries() {
    _lambda = 0.0;
    if (_with_ew_corrections) {
        GluonFusionEWCoefficients ew_coef(_model);
        _lambda = abs(ew_coef.LO());
    }
    // note: here we set the coeff. of lambda of order a_s^2 to zero.
    // note: this affects the total xs at the 0.6 per mille level
    _ew_lambda_series = AsSeries(1,
                                 _lambda,
                                 _lambda * 7./6.,
                                 _lambda * 0.0);
}



void InputParameters::CreateInputInformation() {
    stringstream info;
    info << "-----------------------" << endl;
    info << "Input used in this run:" << endl;
    info << "a_s(m_z) = " << _as_at_mz << endl;
    info << "a_s(mu_r) = " << _as_over_pi*consts::Pi << endl;
    info << "mh = " << _model.higgs.m()
    << "\t (" << _model.higgs.scheme() << ")" << endl;
    info << "mt = " << _model.top.m()
    << "\t (" << _model.top.scheme() << ")" << endl;
    info << "mb = " << _model.bottom.m()
    << "\t (" << _model.bottom.scheme() << ")" << endl;
    info << "mc = " << _model.charm.m()
    << "\t (" << _model.charm.scheme() << ")" << endl;
    info << "ew factor (2*lambda) = " << 2. * _lambda << endl;
    _input_information = info.str();
}


int InputParameters::DetermineQCDPerturbativeOrder(const string& order) {
    if (order == "LO")  return 2;
    if (order == "NLO") return 3;
    if (order == "NNLO") return 4;
    if (order == "N3LO") return 5;
    cout << "\nFailed to determine int(perturbative order) for input: "
    << order << endl;
    exit(EXIT_FAILURE);
}



void InputParametersForThresRes::Configure(const InputParameters & IP){
    _wc = IP._wc;
    _prefactor = IP._prefactor;
    _tau = IP._tau;
    _mur = IP._mur;
    _muf = IP._muf;
    _as_over_pi = IP._as_over_pi;
    _log_muf_over_mh_sq = IP._log_muf_over_mh_sq;
    _log_mur_over_muf_sq = IP._log_mur_over_muf_sq;

}
