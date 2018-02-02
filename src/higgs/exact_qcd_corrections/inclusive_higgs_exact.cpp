
#include "inclusive_higgs_exact.h"


#include "higgs_eft.h"



void InclusiveHiggsExact::SetUpContributions(){
    Channel* gg_channel=new Channel("gg");
    Channel* qg_channel=new Channel("qg");
    Channel* qqbar_channel=new Channel("qqbar");
    AsSeries a_s(1,1.);
    AsSeries a_s_sq(2,1.);
    AsSeries a_s_cube(3,1.);
    
    // computing the exact LO
    //GluonFusionExactCoefficients exact(_input._model);
    
    h_exact::ExactRegularCalculator* exact = new h_exact::ExactRegularCalculator(&(_input._model));
    exact->compute_soft();
    
    
    gg_channel->AddTerm("LO delta exact",
                        a_s_sq * exact->born_squared(),
                        new HiggsExact_gg_delta);
    
    if (_input._int_qcd_perturbative_order>2)
    {
        // the NLO delta exact matrix element in gg channel
        
        gg_channel->AddTerm("NLO delta exact",
                            a_s_cube*exact->soft(),
                            new HiggsExact_gg_delta);
        
        // the NLO D0 exact matrix element in gg channel (which is the same as in the effective but multiplied by the exact born).
        // Before having a heart attack here, remember that n_D0_* returns
        //  an AsSeries that starts at a_s^1, so the total a_s power
        //  here is 3 (as it should)
        
        
        gg_channel->AddTerm("NLO D0 exact",a_s_sq*exact->born_squared()
                                *( HEFT::n_D0_log_muf(_input._log_muf_over_mh_sq) + HEFT::n_D0_at_mh()),
                                new HiggsExact_gg_plus(0)
                                );
        
        
        
        gg_channel->AddTerm("NLO D1 exact",a_s_sq*exact->born_squared()
                                *( HEFT::n_D1_log_muf(_input._log_muf_over_mh_sq) + HEFT::n_D1_at_mh()),
                                new HiggsExact_gg_plus(1)
                                );
        
        // the exact real gluon gluon
        
        gg_channel->AddTerm("NLO reg exact",
                            a_s_cube,
                            new HiggsExact_gg_nlo_reg(&(_input._model))
                            );
        
        // the exact reg quark gluon
        
        qg_channel->AddTerm("NLO reg exact",
                            a_s_cube,
                            new HiggsExact_qg_nlo_reg(&(_input._model))
                            
                            );
        // the exact reg q qbar
        
        qqbar_channel->AddTerm("NLO reg exact",
                            a_s_cube,
                            new HiggsExact_qqbar_nlo_reg(&(_input._model))
                            );
        
    }
    
    _channels.push_back(gg_channel);
    _channels.push_back(qg_channel);
    _channels.push_back(qqbar_channel);
    
    delete exact;
    
}

void InclusiveHiggsExact::EvaluatePostVegasDependingOnProcess(){
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
            //truncate to the requested power: here it must be maximum 3
            if (_input._int_qcd_perturbative_order>2)
                _qcd_result.Truncate(3);
            else
                _qcd_result.Truncate(2);
            
            the_term->StoreQCDResult(_qcd_result);
            
        }
    }
}

HiggsExact_gg_delta::HiggsExact_gg_delta(){
    set_dimensions(2); // for cuhre
    _channel = "gg";
    _name="gg_delta_exact";
   
}

HiggsExact_gg_plus::HiggsExact_gg_plus(int m){
    set_dimensions(2);
    _channel = "gg";
    _name="gg_D"+to_string(m)+"_exact";

    _log_power=m;
    //set_integration_method("vegas");
}

HiggsExact_gg_nlo_reg::HiggsExact_gg_nlo_reg(CModel* model){
    _model=model;
    set_dimensions(3);
    _channel = "gg";
    _name="gg_reg_nlo_exact";
    
    set_rel_accuracy_multiplier(5.);
    _hard_only = false ;
    _calculator = new h_exact::ExactRegularCalculator(model);
    //set_integration_method("vegas");
}

void HiggsExact_gg_nlo_reg::compute_matrix_element_Lf_coefficients(const double* xx){
    if (not(_hard_only)){
        _me_lf_coeffs.clear();
        _me_lf_coeffs.push_back(_calculator->gg_reg_Log0(xx[1],xx[2]));
        _me_lf_coeffs.push_back(_calculator->gg_reg_Log1(xx[1],xx[2]));
       // _me_lf_coeffs.push_back(h_exact::gg_reg_exact_nlo_L0(xx[1],xx[2],*_model));
       // _me_lf_coeffs.push_back(h_exact::gg_reg_exact_nlo_L1(xx[1],xx[2],*_model));
    }
    //: hard only
    else{
        _me_lf_coeffs.clear();
        _me_lf_coeffs.push_back(_calculator->gg_reg_Log0_hard(xx[1],xx[2]));
        // the entire Log[muf/mh] is considered soft-collinear
        _me_lf_coeffs.push_back(0.0);
    }
}


HiggsExact_qg_nlo_reg::HiggsExact_qg_nlo_reg(CModel* model){
    _model=model;
    set_dimensions(3);
    _channel = "qg";
    _name="qg_reg_nlo_exact";
    
    set_rel_accuracy_multiplier(10.);
    _hard_only = false;
    _calculator = new h_exact::ExactRegularCalculator(model);

    //set_integration_method("vegas");
}

void HiggsExact_qg_nlo_reg::compute_matrix_element_Lf_coefficients(const double* xx){
    if ( not(_hard_only) ){
        _me_lf_coeffs.clear();
        _me_lf_coeffs.push_back(_calculator->qg_reg_Log0(xx[1],xx[2]));
        _me_lf_coeffs.push_back(_calculator->qg_reg_Log1(xx[1],xx[2]));
    }
    //: hard only
    else{
        _me_lf_coeffs.clear();
        _me_lf_coeffs.push_back(_calculator->qg_reg_Log0_hard(xx[1],xx[2]));
        // the entire Log[muf/mh] is considered soft-collinear
        _me_lf_coeffs.push_back(0.0);
    }
}


HiggsExact_qqbar_nlo_reg::HiggsExact_qqbar_nlo_reg(CModel* model){
    _model=model;
    set_dimensions(3);
    _channel = "qqbar";
    _name="qqbar_reg_nlo_exact";
    
    set_rel_accuracy_multiplier(10.);
    _calculator = new h_exact::ExactRegularCalculator(model);

    //set_integration_method("vegas");
    
}

void HiggsExact_qqbar_nlo_reg::compute_matrix_element_Lf_coefficients(const double* xx){
    _me_lf_coeffs.clear();
    _me_lf_coeffs.push_back(_calculator->qqb_reg_Log0(xx[1],xx[2]));
}



