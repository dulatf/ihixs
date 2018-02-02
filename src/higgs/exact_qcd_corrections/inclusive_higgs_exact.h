#ifndef INCLUSIVE_HIGGS_EXACT_H
#define INCLUSIVE_HIGGS_EXACT_H

#include "inclusive_process.h"
#include "nlo_exact_matrix_elements.h"//: for h_exact namespace


class InclusiveHiggsExact: public InclusiveProcess{
public:
    InclusiveHiggsExact(const UserInterface& UI):InclusiveProcess(UI){_name = "Higgs Production in gluon fusion: Exact (quark mass effects retained)";}
    void SetUpContributions();
    void EvaluatePostVegasDependingOnProcess();
private:
};


class HiggsExact_gg_delta: public LuminosityIntegralAbsolutelyGenericDelta
{
public:
    HiggsExact_gg_delta();
};

class HiggsExact_gg_plus: public LuminosityIntegralAbsolutelyGenericPlus
{
public:
    HiggsExact_gg_plus(int m);
};

class HiggsExact_gg_nlo_reg: public LuminosityIntegralAbsolutelyGenericReg
{
public:
    HiggsExact_gg_nlo_reg(CModel* model);
    void compute_matrix_element_Lf_coefficients(const double* xx);
    void SetHardOnlyFlag(bool flag){_hard_only = flag;}
    
private:
    CModel* _model;
    h_exact::ExactRegularCalculator* _calculator;
    bool _hard_only;


    
};


class HiggsExact_qg_nlo_reg: public LuminosityIntegralAbsolutelyGenericReg
{
public:
    HiggsExact_qg_nlo_reg(CModel* model);
    void compute_matrix_element_Lf_coefficients(const double* xx);
    void SetHardOnlyFlag(bool flag){_hard_only = flag;}

private:
    CModel* _model;
    h_exact::ExactRegularCalculator* _calculator;

    bool _hard_only;
};

class HiggsExact_qqbar_nlo_reg: public LuminosityIntegralAbsolutelyGenericReg
{
public:
    HiggsExact_qqbar_nlo_reg(CModel* model);
    h_exact::ExactRegularCalculator* _calculator;
    void compute_matrix_element_Lf_coefficients(const double* xx);
private:
    CModel* _model;
    
};






#endif
