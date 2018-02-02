
#include "inclusive_higgs_eft.h"


void InclusiveHiggsEFT::SetUpContributions(){
    
    
    
    AddGG();
    

    AddQG();
    AddQQBAR();
    AddQQ();
    AddQ1Q2();
    
    
}


void InclusiveHiggsEFT::AddGG(){
    AsSeries a_s(1,1.);
    AsSeries a_s_sq(2,1.);
    AsSeries a_s_cube(3,1.);
    
    
    Channel* gg_channel=new Channel("gg");
    
    double _log_muf_mh_sq = _input._log_muf_over_mh_sq;
    // constructing terms
    
    gg_channel->AddTerm("delta ",
     HEFT::n_delta_at_mh()
     +HEFT::n_delta_log_muf(_log_muf_mh_sq),
     new HiggsEFTGGDelta);
    
    if (_input._int_qcd_perturbative_order>2){
        
        gg_channel->AddTerm("D0 ",
         HEFT::n_D0_at_mh()
         +HEFT::n_D0_log_muf(_log_muf_mh_sq),
         new HiggsEFTGGPlus(0));
         gg_channel->AddTerm("D1 ",
         HEFT::n_D1_at_mh()
         +HEFT::n_D1_log_muf(_log_muf_mh_sq),
         new HiggsEFTGGPlus(1));
        
        
        
        
        
        gg_channel->AddTerm("NLO reg ",a_s,new HiggsEFTGGNLOReg);
    }
    
    if (_input._int_qcd_perturbative_order>3){
        gg_channel->AddTerm("D2 ",
                            HEFT::n_D2_at_mh()
                            +HEFT::n_D2_log_muf(_log_muf_mh_sq),
                            new HiggsEFTGGPlus(2));
        gg_channel->AddTerm("D3 ",
                            HEFT::n_D3_at_mh()
                            +HEFT::n_D3_log_muf(_log_muf_mh_sq),
                            new HiggsEFTGGPlus(3));
        
        gg_channel->AddTerm("NNLO reg ",a_s_sq,new HiggsEFTGGNNLOReg);
        
    }
    if (_input._int_qcd_perturbative_order>4){
        gg_channel->AddTerm("D4 ",
                            HEFT::n_D4_at_mh()
                            +HEFT::n_D4_log_muf(_log_muf_mh_sq),
                            new HiggsEFTGGPlus(4));
        gg_channel->AddTerm("D5 ",
                            HEFT::n_D5_at_mh(),
                            new HiggsEFTGGPlus(5));
        gg_channel->AddTerm("N3LO reg ",a_s_cube,new HiggsEFTGGN3LOReg);
        
    }
    
    _channels.push_back(gg_channel);
    
}

void InclusiveHiggsEFT::AddQG(){
    Channel* qg_channel=new Channel("qg");
    
    AsSeries a_s(1,1.);
    AsSeries a_s_sq(2,1.);
    AsSeries a_s_cube(3,1.);
    
    if (_input._int_qcd_perturbative_order>2){
        qg_channel->AddTerm("NLO reg ",a_s,new HiggsEFTQGNLOReg);
        
    }
    
    if (_input._int_qcd_perturbative_order>3){
        qg_channel->AddTerm("NNLO reg ",a_s_sq,new HiggsEFTQGNNLOReg);
        
    }
    
    if (_input._int_qcd_perturbative_order>4){
        qg_channel->AddTerm("N3LO reg ",a_s_cube,new HiggsEFTQGN3LOReg);
    }
    _channels.push_back(qg_channel);
    
}

void InclusiveHiggsEFT::AddQQBAR(){
    Channel* qqbar_channel=new Channel("qqbar");
    
    AsSeries a_s(1,1.);
    AsSeries a_s_sq(2,1.);
    AsSeries a_s_cube(3,1.);
    if (_input._int_qcd_perturbative_order>2){
        qqbar_channel->AddTerm("NLO reg ",a_s,new HiggsEFTQQBARNLOReg);
    }
    
    if (_input._int_qcd_perturbative_order>3){
        qqbar_channel->AddTerm("NNLO reg ",a_s_sq,new HiggsEFTQQBARNNLOReg);
        
    }
    
    if (_input._int_qcd_perturbative_order>4){
        
        qqbar_channel->AddTerm("N3LO reg ",a_s_cube,new HiggsEFTQQBARN3LOReg);
        
    }
    _channels.push_back(qqbar_channel);
    
}

void InclusiveHiggsEFT::AddQQ(){
    Channel* qq_channel=new Channel("qq");
    
    AsSeries a_s(1,1.);
    AsSeries a_s_sq(2,1.);
    AsSeries a_s_cube(3,1.);
    if (_input._int_qcd_perturbative_order>3){
        qq_channel->AddTerm("NNLO reg ",a_s_sq,new HiggsEFTQQNNLOReg);
        
    }
    
    if (_input._int_qcd_perturbative_order>4){
        qq_channel->AddTerm("N3LO reg ",a_s_cube,new HiggsEFTQQN3LOReg);
        
    }
    _channels.push_back(qq_channel);
    
}

void InclusiveHiggsEFT::AddQ1Q2(){
    Channel* q1q2_channel=new Channel("q1q2");
    
    AsSeries a_s(1,1.);
    AsSeries a_s_sq(2,1.);
    AsSeries a_s_cube(3,1.);
    if (_input._int_qcd_perturbative_order>3){
        q1q2_channel->AddTerm("NNLO reg ",a_s_sq,new HiggsEFTQ1Q2NNLOReg);
        
    }
    
    if (_input._int_qcd_perturbative_order>4){
         q1q2_channel->AddTerm("N3LO reg ",a_s_cube,new HiggsEFTQ1Q2N3LOReg);
        
    }
    _channels.push_back(q1q2_channel);
}


void InclusiveHiggsEFT::EvaluatePostVegasDependingOnProcess(){
    for (int i=0;i<_channels.size();i++){
        for (int j=0;j<_channels[i]->NumberOfTerms();j++){
            SigmaTerm* the_term = _channels[i]->GiveTermPtr(j);
            
            AsSeries res = the_term->PostVegasResult();
            
            // multiplying by universal prefactor (equal to effective theory born divided by a_s^2)
            res = res * _input._prefactor;
            
            // multiply by C^2
            AsSeries _qcd_result = res * _input._wc.c() * _input._wc.c();
            // truncate to order a_s^5
            _qcd_result.Truncate(5);
            
            
            //: electroweak corrections
            // set up the weak wilson coeff: (C+lambda*(1+a_s 7/6 + a_s^2 10)^2 -C^2)
            // and multiply with it to get the ew contribution
            AsSeries _ew_result = res *  (2. * _input._wc.c() * _input._ew_lambda_series
                                          + _input._ew_lambda_series * _input._ew_lambda_series);
            // truncate ew to otder a_s^5
            _ew_result.Truncate(5);
            
            if (abs(_input._log_mur_over_muf_sq) > 1e-15)
                // mu_r evolution
            {
                _qcd_result = MurEvolution(_qcd_result,-_input._log_mur_over_muf_sq);
                _ew_result = MurEvolution(_ew_result,-_input._log_mur_over_muf_sq);
            }
            
            
            
            // multiply each term with a_s/pi to the proper power
            _qcd_result.MultiplyAs(_input._as_over_pi);
            _ew_result.MultiplyAs(_input._as_over_pi);
            //truncate to the requested power
            _qcd_result.Truncate(_input._int_qcd_perturbative_order);
            _ew_result.Truncate(_input._int_qcd_perturbative_order);
            
            the_term->StoreQCDResult(_qcd_result);
            the_term->StoreEWResult(_ew_result);
            
            
            // WC^2*n implementation (beyond fixed order)
            
            AsSeries eta = res;
            AsSeries wc = _input._wc.c();
            
            eta.ShiftStartingExponentBy(2);
            wc.ShiftStartingExponentBy(-1);
            
            if (abs(_input._log_mur_over_muf_sq) > 1e-15)
                // mu_r evolution
            {
                eta = GenericMurEvolution(eta,-_input._log_mur_over_muf_sq);
                wc = GenericMurEvolution(wc,-_input._log_mur_over_muf_sq);
            }
            
            
            eta.MultiplyAs(_input._as_over_pi);
            wc.MultiplyAs(_input._as_over_pi);
            
            eta.Truncate(5);
            wc.Truncate(3);
            
            the_term->StoreEta(eta);
            the_term->StoreWC(wc);
            
            
        }
    }
}

AsSeries InclusiveHiggsEFT::EwCorrections()
{
    AsSeries res(2,0.0);
    for (int j=0;j<_channels.size();j++)
    {
        res=res+_channels[j]->EWResult();
    }
    
    return res;
}

AsSeries InclusiveHiggsEFT::Eta(){
    AsSeries res;
    for (int i=0;i<_channels.size();i++){
        for (int j=0;j<_channels[i]->NumberOfTerms();j++){
            SigmaTerm* the_term = _channels[i]->GiveTermPtr(j);
            res = res + the_term->Eta();
        }
    }
    return res;
}

ResultPair InclusiveHiggsEFT::R_LO_x_EFT(int order) {
    if (order < 2) {
        cout << "Error in InclusiveHiggsEFT::R_LO_x_EFT(int order): we count orders from LO=2 and you asked for " << order << endl;
        exit(1);
    }
    if (order > 5) {
        cout << "Error in InclusiveHiggsEFT::R_LO_x_EFT(int order): we count orders until N3LO=5 and you asked for " << order << " which we haven't calculated yet." << endl;
        exit(1);
    }
    
    AsSeries eft_res = QCDCorrections();
    ResultPair res;
    for (int i=2; i < order+1; i++) {
        res = res + eft_res.term_of_order(i);
    }
    return RescalingCoefficient() * res ;
}






