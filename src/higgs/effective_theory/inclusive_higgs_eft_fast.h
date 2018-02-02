#ifndef INCLUSIVE_HIGGS_EFT_FAST_H
#define INCLUSIVE_HIGGS_EFT_FAST_H

#include "inclusive_process.h"


class InclusiveHiggsEFTFast: public InclusiveProcess{
public:
    InclusiveHiggsEFTFast(const UserInterface& UI):InclusiveProcess(UI){_name = "Higgs Production in gluon fusion: EFT fast version";}
    void SetUpContributions();
    void EvaluatePostVegasDependingOnProcess();

    
    double RescalingCoefficient(){return _input._top_only_LO_coefficient;}
    // EW corrections as a series in a_s
    AsSeries EwCorrections();
    
    AsSeries WC(){return _channels[0]->GiveTermPtr(0)->WC();}
    AsSeries Eta();
    
    vector<double> ResultVector();
    vector<double> ErrorVector();
    
};

#include "higgs_eft.h"



class HiggsEFTGGfast: public LuminosityIntegral{
public:
    HiggsEFTGGfast();
    vector<double> evaluateIntegral(const double* xx);
    void SetInput(InputParameters* inputptr){_input=inputptr;}

    void EvaluateDeltaAndPlusCoeffs();
private:
    // all input for the xs
    InputParameters* _input;
    AsSeries _delta;
    AsSeries _d0;
    AsSeries _d1;
    AsSeries _d2;
    AsSeries _d3;
    AsSeries _d4;
    AsSeries _d5;
    bool _delta_and_plus_evaluated;
    double _K;
    double _exact_LO_delta;
    double _exact_NLO_delta;
    double _exact_d0;
    double _exact_d1;
    CModel* _model_ptr;
    HEFT::N3LORegEvaluator evaluator;
};

class HiggsEFTQGfast: public LuminosityIntegral{
public:
    HiggsEFTQGfast();
    vector<double> evaluateIntegral(const double* xx);
    void SetInput(InputParameters* inputptr){_input=inputptr;}
    void SetK();
private:
    // all input for the xs
    InputParameters* _input;
    double _K;
    bool _k_computed;
    CModel* _model_ptr;
    HEFT::N3LORegEvaluator evaluator;

};

class HiggsEFTQQBARfast: public LuminosityIntegral{
public:
    HiggsEFTQQBARfast();
    vector<double> evaluateIntegral(const double* xx);
    void SetInput(InputParameters* inputptr){_input=inputptr;}
    void SetK();
private:
    // all input for the xs
    InputParameters* _input;
    double _K;
    bool _k_computed;
    CModel* _model_ptr;
    HEFT::N3LORegEvaluator evaluator;

};


class HiggsEFTQQfast: public LuminosityIntegral{
public:
    HiggsEFTQQfast();
    vector<double> evaluateIntegral(const double* xx);
    void SetInput(InputParameters* inputptr){_input=inputptr;}

    void SetK();

private:
    // all input for the xs
    InputParameters* _input;
    double _K;
    bool _k_computed;
    HEFT::N3LORegEvaluator evaluator;

};

class HiggsEFTQ1Q2fast: public LuminosityIntegral{
public:
    HiggsEFTQ1Q2fast();
    vector<double> evaluateIntegral(const double* xx);
    void SetInput(InputParameters* inputptr){_input=inputptr;}
    
    void SetK();
private:
    // all input for the xs
    InputParameters* _input;
    double _K;
    bool _k_computed;
    HEFT::N3LORegEvaluator evaluator;

};

#endif
