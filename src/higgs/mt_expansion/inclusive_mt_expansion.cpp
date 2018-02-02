#include "inclusive_mt_expansion.h"
#include "nlo_exact_matrix_elements.h"
#include "higgs_eft.h"
#include "mt_expansion_gg.h"
#include "mt_expansion_qg.h"

//------------------ InclusiveHiggsMtExpansion

void InclusiveHiggsMtExpansion::SetUpContributions(){

    Channel* gg_channel=new Channel("gg");
    Channel* qg_channel=new Channel("qg");
    //Channel* qqbar_channel=new Channel("qqbar");
    
    
    AsSeries a_s_sq(2,1.);
    AsSeries a_s_cube(3,1.);
    AsSeries a_s_fourth(4,1.);
    
    
    
    const double rho =  pow(_input._model.higgs.m()/_input._model.top.m(),2.);
    
    // computing the exact LO top only
    _input._model.bottom.set_Y(0.0);
    _input._model.charm.set_Y(0.0);
    
    
    
    h_exact::ExactRegularCalculator exact(&(_input._model));
    
    
    // the exact LO delta in gg channel
    // the NLO delta exact matrix element in gg channel
    //gg_channel->AddDeltaTerm("LO delta exact",a_s_sq*exact.born_squared());
    

    if (_input._int_qcd_perturbative_order>3)
    {
        
        
        gg_channel->AddTerm("NNLO 1/mt : delta",
                            a_s_fourth
                            *(MTEXP::gg_nnlo_delta(rho,
                                                  _input._log_muf_over_mh_sq,
                                                  _input._log_muf_over_mt_sq
                                                   )
//                                -
//                              MTEXP::gg_nnlo_delta(0.0,
//                                                   _input._log_muf_over_mh_sq,
//                                                   _input._log_muf_over_mt_sq
//                                                   )
                                )
                            *exact.born_squared(),
                            new HiggsMtExpansion_gg_nnlo_delta
                            );
        
        gg_channel->AddTerm("NNLO 1/mt: D0",
                            a_s_fourth*(MTEXP::gg_nnlo_D0(rho,_input._log_muf_over_mh_sq)
                            //-MTEXP::gg_nnlo_D0(0.0,_input._log_muf_over_mh_sq)
                            )*exact.born_squared(),
                            new HiggsMtExpansion_gg_nnlo_plus(0));
        
        gg_channel->AddTerm("NNLO 1/mt: D1",
                            a_s_fourth*(MTEXP::gg_nnlo_D1(rho,_input._log_muf_over_mh_sq)
                            //-MTEXP::gg_nnlo_D1(0.0,_input._log_muf_over_mh_sq)
                            )*exact.born_squared(),
                            new HiggsMtExpansion_gg_nnlo_plus(1));
        
        gg_channel->AddTerm("NNLO 1/mt: D2",
                            a_s_fourth*(MTEXP::gg_nnlo_D2(rho,_input._log_muf_over_mh_sq)
                            //-MTEXP::gg_nnlo_D2(0.0,_input._log_muf_over_mh_sq)
                            )*exact.born_squared(),
                            new HiggsMtExpansion_gg_nnlo_plus(2));
        
        gg_channel->AddTerm("NNLO 1/mt: D3",
                            a_s_fourth*(MTEXP::gg_nnlo_D3(rho,_input._log_muf_over_mh_sq)
                            //-MTEXP::gg_nnlo_D3(0.0,_input._log_muf_over_mh_sq)
                            )*exact.born_squared(),
                            new HiggsMtExpansion_gg_nnlo_plus(3));
        
        gg_channel->AddTerm("NNLO 1/mt: regular",
                            a_s_fourth*exact.born_squared(),
                            new HiggsMtExpansion_gg_nnlo_reg(rho));
//        gg_channel->AddTerm("NNLO 1/mt: -eft regular ",
//                            a_s_fourth*(-1.)*exact.born_squared(),
//                            new HiggsMtExpansion_gg_nnlo_reg(0.0));
        
        
        //qg channel
        qg_channel->AddTerm("NNLO 1/mt: regular",
                            a_s_fourth*exact.born_squared(),
                            new HiggsMtExpansion_qg_nnlo_reg(rho));
//        qg_channel->AddTerm("NNLO 1/mt: -eft regular",
//                            a_s_fourth*(-1.)*exact.born_squared(),
//                            new HiggsMtExpansion_qg_nnlo_reg(0.0));
        
        
    }
    
    
    _channels.push_back(gg_channel);
    _channels.push_back(qg_channel);
    
    
}

void InclusiveHiggsMtExpansion::EvaluatePostVegasDependingOnProcess(){
    for (int i=0;i<_channels.size();i++){
        for (int j=0;j<_channels[i]->NumberOfTerms();j++){
            SigmaTerm* the_term = _channels[i]->GiveTermPtr(j);
            
            AsSeries res = the_term->PostVegasResult();
            
            // multiplying by universal prefactor (equal to effective theory born)
            AsSeries _qcd_result = res * _input._prefactor;
            
            if (abs(_input._log_mur_over_muf_sq) > 1e-15)
                // mu_r evolution
            {
                _qcd_result = MurEvolution(_qcd_result,-_input._log_mur_over_muf_sq);
            }
            
            
            
            // multiply each term with a_s/pi to the proper power
            _qcd_result.MultiplyAs(_input._as_over_pi);
            //truncate to the requested power
            _qcd_result.Truncate(_input._int_qcd_perturbative_order);
            
            the_term->StoreQCDResult(_qcd_result);
            
        }
    }
}

ResultPair InclusiveHiggsMtExpansion::delta_mt_gluon_gluon(const ResultPair gg_eft_rescaled_nnlo){
    ResultPair gg_mt = GiveChannel("gg")->QCDResult().term_of_order(4);
    
    return gg_mt - gg_eft_rescaled_nnlo;
}


ResultPair InclusiveHiggsMtExpansion::delta_mt_quark_gluon(const ResultPair qg_eft_rescaled_nnlo){
    ResultPair qg_mt = GiveChannel("qg")->QCDResult().term_of_order(4);
    
    return qg_mt - qg_eft_rescaled_nnlo;
}

