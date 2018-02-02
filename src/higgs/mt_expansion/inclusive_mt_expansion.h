#ifndef INCLUSIVE_MT_EXPANSION
#define INCLUSIVE_MT_EXPANSION


#include "inclusive_process.h"



class InclusiveHiggsMtExpansion : public InclusiveProcess{
public:
    InclusiveHiggsMtExpansion(const UserInterface& UI):InclusiveProcess(UI){
            _name = "Higgs Production in gluon fusion: m_top expansion";}
    void SetUpContributions();
    void EvaluatePostVegasDependingOnProcess();
    ResultPair delta_mt_gluon_gluon(const ResultPair gg_eft_rescaled_nnlo);
    ResultPair delta_mt_quark_gluon(const ResultPair qg_eft_rescaled_nnlo);
};




#endif