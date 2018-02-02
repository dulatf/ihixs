
#include "inclusive_higgs_eft_fast.h"
#include "constants.h"
#include "nlo_exact_matrix_elements.h"//: for h_exact namespace
//: we need this because we allow the fast mode to run with nlo_qcd_exact effects



AsSeries SpecialMurEvolution(const AsSeries& expansion, const double& L){
    const double b0 = consts::beta_zero;
    const double b1 = consts::beta_one;
    const double b2 = consts::beta_two;
    
    const ResultPair c2=expansion.term_of_order(2);
    const ResultPair c3=expansion.term_of_order(3);
    const ResultPair c4=expansion.term_of_order(4);
    const ResultPair c5=expansion.term_of_order(5);
    
    AsSeries res2(2,c2);
    AsSeries res3(3,c3 - 2* b0* c2* L);
    AsSeries res4(4,c4 + (-2*b1*c2-3*b0*c3)*L + 3*c2*pow(b0*L,2));
    AsSeries res5(5,c5 + (-2 * b2 * c2 - 3 * b1 * c3 - 4 * b0 * c4) * L
                  + (7 * b0 * b1*  c2 + 6 * pow(b0,2.) *c3)*pow(L,2.)
                  + (-4 *pow(b0,3.)* c2)*pow(L,3.));
    
    return res2+res3+res4+res5;
}


void InclusiveHiggsEFTFast::SetUpContributions(){
    
    
    AsSeries one(0,1.);
    Channel* gg_channel=new Channel("gg");
    HiggsEFTGGfast* gg = new HiggsEFTGGfast;
    gg->SetInput(&_input);
    gg_channel->AddTerm("total",one,gg);
    
    
    Channel* qg_channel=new Channel("qg");
    HiggsEFTQGfast* qg = new HiggsEFTQGfast;
    qg->SetInput(&_input);
    qg_channel->AddTerm("total",one,qg);

    Channel* qqbar_channel=new Channel("qqbar");
    HiggsEFTQQBARfast* qqbar = new HiggsEFTQQBARfast;
    qqbar->SetInput(&_input);
    qqbar_channel->AddTerm("total",one,qqbar);
    
    Channel* qq_channel=new Channel("qq");
    HiggsEFTQQfast* qq = new HiggsEFTQQfast;
    qq->SetInput(&_input);
    qq_channel->AddTerm("total",one,qq);
    
    
    Channel* q1q2_channel=new Channel("q1q2");
    HiggsEFTQ1Q2fast* q1q2 = new HiggsEFTQ1Q2fast;
    q1q2->SetInput(&_input);
    q1q2_channel->AddTerm("total",one,q1q2);
    
    _channels.push_back(gg_channel);
    _channels.push_back(qg_channel);
    _channels.push_back(qqbar_channel);
    _channels.push_back(qq_channel);
    _channels.push_back(q1q2_channel);
    
    
}


vector<double> InclusiveHiggsEFTFast::ResultVector(){
    vector<vector<double> > result_vectors_from_channels;
    
    for (int j=0;j<_channels.size();j++){
        result_vectors_from_channels.push_back( _channels[j]->ResultVector());
    }
    
    int number_of_components=result_vectors_from_channels[0].size();
    
    vector<double> res(number_of_components,0.0);
    for (int i=0;i<number_of_components;i++){
        for (int j=0;j<_channels.size();j++){
            res[i] += result_vectors_from_channels[j][i];
        }
    }
    return res;
}


vector<double> InclusiveHiggsEFTFast::ErrorVector(){
    vector<vector<double> > result_vectors_from_channels;
    
    for (int j=0;j<_channels.size();j++){
        result_vectors_from_channels.push_back( _channels[j]->ErrorVector());
    }
    
    int number_of_components=result_vectors_from_channels[0].size();
    
    vector<double> res(number_of_components,0.0);
    for (int i=0;i<number_of_components;i++){
        for (int j=0;j<_channels.size();j++){
            res[i] += pow(result_vectors_from_channels[j][i],2.);
        }
        res[i] = sqrt(res[i]);
    }
    return res;
}


void InclusiveHiggsEFTFast::EvaluatePostVegasDependingOnProcess(){
    for (int i=0;i<_channels.size();i++){
        for (int j=0;j<_channels[i]->NumberOfTerms();j++){
            // here we just store the result of each term inside each term!
            // Everything has already been computed, so there is no post-vegas evaluation
            // the only reason for this to exist is compliance with the detailed mode
            // structure, where every term at every order is a different vegas run
            SigmaTerm* the_term = _channels[i]->GiveTermPtr(j);
            the_term->StoreQCDResult(the_term->PostVegasResult());
        }
    }
}

AsSeries InclusiveHiggsEFTFast::EwCorrections()
{
    AsSeries res(2,0.0);
    for (int j=0;j<_channels.size();j++)
    {
        res=res+_channels[j]->EWResult();
    }
    
    return res;
}

AsSeries InclusiveHiggsEFTFast::Eta(){
    AsSeries res;
    for (int i=0;i<_channels.size();i++){
        for (int j=0;j<_channels[i]->NumberOfTerms();j++){
            SigmaTerm* the_term = _channels[i]->GiveTermPtr(j);
            res = res + the_term->Eta();
        }
    }
    return res;
}




HiggsEFTGGfast::HiggsEFTGGfast(){
    set_dimensions(2);
    _channel = "gg";
    _name="fast";
    _delta_and_plus_evaluated=false;
    
}

void HiggsEFTGGfast::EvaluateDeltaAndPlusCoeffs(){
    //Delta terms
    _delta = HEFT::n_delta_at_mh()
    +HEFT::n_delta_log_muf(_Lf);
    
    //Plus terms
    _d0 = HEFT::n_D0_at_mh() + HEFT::n_D0_log_muf(_Lf);
    _d1 = HEFT::n_D1_at_mh() + HEFT::n_D1_log_muf(_Lf);
    _d2 = HEFT::n_D2_at_mh() + HEFT::n_D2_log_muf(_Lf);
    _d3 = HEFT::n_D3_at_mh() + HEFT::n_D3_log_muf(_Lf);
    _d4 = HEFT::n_D4_at_mh() + HEFT::n_D4_log_muf(_Lf);
    _d5 = HEFT::n_D5_at_mh();
    _delta_and_plus_evaluated = true;
    
    //exact LO delta
    // computing the exact LO
    //GluonFusionExactCoefficients exact(_input->_model);
    h_exact::ExactRegularCalculator exact( &(_input->_model) );
    // the exact LO delta in gg channel
    // the NLO delta exact matrix element in gg channel
    _exact_LO_delta = _input->_prefactor
                    * pow(_input->_as_over_pi,2.)
                    * exact.born_squared();
    _exact_NLO_delta = _input->_prefactor
                    * pow(_input->_as_over_pi,3.)
                    * exact.soft();

    if (abs(_input->_log_mur_over_muf_sq) > 1e-15)
        // mu_r evolution
    {
        _exact_NLO_delta = _exact_NLO_delta - 2.*consts::beta_zero * _exact_LO_delta * (-_input->_log_mur_over_muf_sq);
    }

    _exact_d0=_input->_prefactor
                * pow(_input->_as_over_pi,3.)
                * exact.born_squared()
                        *( HEFT::n_D0_log_muf(_Lf).term_of_order(1).val() + HEFT::n_D0_at_mh().term_of_order(1).val());
    
    
    
    _exact_d1=_input->_prefactor
                * pow(_input->_as_over_pi,3.)
                 *exact.born_squared()
                        *( HEFT::n_D1_log_muf(_Lf).term_of_order(1).val() + HEFT::n_D1_at_mh().term_of_order(1).val());


    _K = _input->_top_only_LO_coefficient;
    _model_ptr = &(_input->_model);
}



vector<double> HiggsEFTGGfast::evaluateIntegral(const double* xx){
    const double x1= xx[0];
    const double z = xx[1];
    double measure = 1./x1;
    
    
    if (not(_delta_and_plus_evaluated)) EvaluateDeltaAndPlusCoeffs();
    AsSeries delta = _delta;
    
    AsSeries plus =   _d0
                    + _d1 * pow(log(1.-z), 1.)
                    + _d2 * pow(log(1.-z), 2.)
                    + _d3 * pow(log(1.-z), 3.)
                    + _d4 * pow(log(1.-z), 4.)
                    + _d5 * pow(log(1.-z), 5.);
    
    // regular terms
    
    //nlo
    const double reg_nlo = HEFT::nlo_reg(z)
                            +HEFT::nlo_reg_L(z) * _Lf;

    //nnlo
    
    // L^0
    const double reg_nnlo =HEFT::gg_n2lo_lzbar0_lz0_no_log(z)
                            +HEFT::gg_n2lo_lzbar0_lz1_no_log(z)
                            +HEFT::gg_n2lo_lzbar0_lz2_no_log(z)
                            +HEFT::gg_n2lo_lzbar0_lz3_no_log(z)
                            +HEFT::gg_n2lo_lzbar1_L_0(z)
                            +HEFT::gg_n2lo_lzbar2_L_0(z)
                            +HEFT::gg_n2lo_lzbar3_L_0(z)
                            +
                            (
                             HEFT::gg_n2lo_lzbar0_lz0_L_1(z)
                            +HEFT::gg_n2lo_lzbar0_lz1_L_1(z)
                            +HEFT::gg_n2lo_lzbar0_lz2_L_1(z)
                            +HEFT::gg_n2lo_lzbar1_L_1(z)
                            +HEFT::gg_n2lo_lzbar2_L_1(z)
                            ) * _Lf
                            +
                            (
                            HEFT::gg_n2lo_lzbar1_L_2(z)
                            +HEFT::gg_n2lo_lzbar0_lz0_L_2(z)
                            +HEFT::gg_n2lo_lzbar0_lz1_L_2(z)
                            ) * pow(_Lf,2.);
    
    const double reg_n3lo = evaluator.n3lo_reg_complete_gg(z,0)
                            -evaluator.n3lo_reg_complete_gg(z,1)* _Lf
                            +evaluator.n3lo_reg_complete_gg(z,2) * pow(_Lf,2.)
                            -evaluator.n3lo_reg_complete_gg(z,3) * pow(_Lf,3.);
    
    
    
    
    AsSeries reg(1,reg_nlo,reg_nnlo,reg_n3lo);
    
    
    
    // multiplying by universal prefactor (equal to effective theory born divided by a_s^2)
    delta = delta * _input->_prefactor;
    plus = plus * _input->_prefactor;
    reg = reg * _input->_prefactor;
    
    // multiply by C^2
    double ew_delta_used;
    double ew_plus_used;
    double ew_reg_used;

       ew_reg_used=0.0;
       ew_plus_used=0.0;
       ew_delta_used=0.0;
    
    delta = delta * _input->_wc.c() * _input->_wc.c();
    plus = plus * _input->_wc.c() * _input->_wc.c();
    reg = reg * _input->_wc.c() * _input->_wc.c();
    
    

    if (abs(_input->_log_mur_over_muf_sq) > 1e-15)
        // mu_r evolution
    {
        delta = SpecialMurEvolution(delta,-_input->_log_mur_over_muf_sq);
        plus = SpecialMurEvolution(plus,-_input->_log_mur_over_muf_sq);
        reg = SpecialMurEvolution(reg,-_input->_log_mur_over_muf_sq);
        
    }
    
    
    
    // multiply each term with a_s/pi to the proper power
    delta.MultiplyAs(_input->_as_over_pi);
    plus.MultiplyAs(_input->_as_over_pi);
    reg.MultiplyAs(_input->_as_over_pi);
    
    //truncate to the requested power
    delta.Truncate(_input->_int_qcd_perturbative_order);
    plus.Truncate(_input->_int_qcd_perturbative_order);
    reg.Truncate(_input->_int_qcd_perturbative_order);
    

    
    double reg_nlo_used;
    double delta_lo_used;
    double delta_nlo_used;

    double plus_nlo_used;
 
        reg_nlo_used = _K*reg.term_of_order(3).val();
        plus_nlo_used =_K*plus.term_of_order(3).val();
        delta_lo_used = _K*delta.term_of_order(2).val();
        delta_nlo_used = _K*delta.term_of_order(3).val();

    //mur<>muf evolution: only affects delta term which is precomputed, up to NLO
    
    double finalDelta =   delta_lo_used
                        + delta_nlo_used
                        + _K*delta.term_of_order(4).val()
                        + _K*delta.term_of_order(5).val()
                        + ew_delta_used;
    
    
    double finalPlus =    plus_nlo_used
                        + _K*plus.term_of_order(4).val()
                        + _K*plus.term_of_order(5).val()
                        + ew_plus_used;
    
    double finalReg =   reg_nlo_used
                      + _K*reg.term_of_order(4).val()
                      + _K*reg.term_of_order(5).val()
                      + ew_reg_used;
    
    vector<double> res;
    vector<double> lumiAt1 = lumi_->luminosity_vector(x1,tau_/x1,_muF);
    vector<double> lumiAtZ =  lumi_->luminosity_vector(x1,tau_/x1/z,_muF);

    for (int i=0;i<_number_of_components;i++){
    res.push_back(  measure * (
                            lumiAt1[i] * finalDelta
                            +(lumiAtZ[i]-lumiAt1[i]) / (1.-z) * finalPlus
                            +lumiAtZ[i] * finalReg
                            )
                        );
    }
    
    

    return res;
}




HiggsEFTQGfast::HiggsEFTQGfast(){
    set_dimensions(2);
    _channel = "qg";
    _name="fast";
    _k_computed=false;
    set_rel_accuracy_multiplier(10.0);

}

void HiggsEFTQGfast::SetK(){
    _K=_input->_top_only_LO_coefficient;
    _k_computed = true;
    _model_ptr = &(_input->_model);

}

vector<double> HiggsEFTQGfast::evaluateIntegral(const double* xx){
    const double x1= xx[0];
    const double z = xx[1];
    //const double lambda = xx[2];
    double measure = 1./x1;
    
    if (not(_k_computed)) SetK();


       // regular terms
    
    //nlo
    const double reg_nlo =    HEFT::qg_nlo_r_L0(xx[1])
                            + HEFT::qg_nlo_r_L1(xx[1])* _Lf;
    
    //nnlo
    
    const double reg_nnlo =HEFT::qg_n2lo_r_lz0_L0(z)
    +HEFT::qg_n2lo_r_lz1_L0(z)
    +HEFT::qg_n2lo_r_lz2_L0(z)
    +HEFT::qg_n2lo_r_lz3_L0(z)
    +
    (
     HEFT::qg_n2lo_r_lz0_L1(z)
     +HEFT::qg_n2lo_r_lz1_L1(z)
     +HEFT::qg_n2lo_r_lz2_L1(z)
     ) * _Lf
    +
    (
     HEFT::qg_n2lo_r_lz0_L2(z)
     +HEFT::qg_n2lo_r_lz1_L2(z)
     ) * pow(_Lf,2.);
    
    const double reg_n3lo = evaluator.n3lo_reg_complete_qg(z,0)
    -evaluator.n3lo_reg_complete_qg(z,1) * _Lf
    +evaluator.n3lo_reg_complete_qg(z,2) * pow(_Lf,2.)
    -evaluator.n3lo_reg_complete_qg(z,3) * pow(_Lf,3.);

    
    AsSeries reg(1,reg_nlo,reg_nnlo,reg_n3lo);
    
    // multiplying by universal prefactor (equal to effective theory born divided by a_s^2)
   
    reg=reg * _input->_prefactor;
    
    // multiply by C^2
    double reg_ew_used=0.0;

        reg_ew_used = 0.0;
    
    reg = reg * _input->_wc.c() * _input->_wc.c();
    
    
    if (abs(_input->_log_mur_over_muf_sq) > 1e-15)
    {
    
        reg = SpecialMurEvolution(reg,-_input->_log_mur_over_muf_sq);
        
    }
    
    // multiply each term with a_s/pi to the proper power
    
    reg.MultiplyAs(_input->_as_over_pi);
    
    //truncate to the requested power
    
    reg.Truncate(_input->_int_qcd_perturbative_order);
    
    
    
    //reg
    double reg_nlo_used;

        reg_nlo_used = _K*reg.term_of_order(3).val();
    
    double finalReg = reg_nlo_used
                     + _K*reg.term_of_order(4).val()
                     + _K*reg.term_of_order(5).val()
                     + reg_ew_used;
    
    
    vector<double> res;
    vector<double> lumiAtZ =  lumi_->luminosity_vector(x1,tau_/x1/z,_muF);
    
    for (int i=0;i<_number_of_components;i++){
        res.push_back(  measure * lumiAtZ[i] * finalReg);
    }
    
    
    return res;
}


HiggsEFTQQBARfast::HiggsEFTQQBARfast(){
    set_dimensions(2);
    _channel = "qqbar";
    _name="fast";
    _k_computed=false;
    set_rel_accuracy_multiplier(10.0);
}

void HiggsEFTQQBARfast::SetK(){
    _K=_input->_top_only_LO_coefficient;
    _k_computed = true;
    _model_ptr = &(_input->_model);
}

vector<double> HiggsEFTQQBARfast::evaluateIntegral(const double* xx){
    const double x1= xx[0];
    const double z = xx[1];
//    const double lambda = xx[2];
    double measure = 1./x1;
    
    if (not(_k_computed)) SetK();

    
    // regular terms
    
    //nlo
    const double reg_nlo =    HEFT::qqbar_nlo_r_L0(z);
    
    //nnlo
    
    // L^0
    const double reg_nnlo =HEFT::qqb_nnlo_r_lz0_L0(z)
    +HEFT::qqb_nnlo_r_lz1_L0(z)
    +HEFT::qqb_nnlo_r_lz2_L0(z)
    +
    (
     HEFT::qqb_nnlo_r_lz0_L1(z)
     +HEFT::qqb_nnlo_r_lz1_L1(z)
     ) * _Lf
    +
    (
     HEFT::qqb_nnlo_r_lz0_L2(z)
     ) * pow(_Lf,2.);
    
    const double reg_n3lo = evaluator.n3lo_reg_complete_qqbar(z,0)
    -evaluator.n3lo_reg_complete_qqbar(z,1) * _Lf
    +evaluator.n3lo_reg_complete_qqbar(z,2) * pow(_Lf,2.)
    -evaluator.n3lo_reg_complete_qqbar(z,3) * pow(_Lf,3.);

    
    AsSeries reg(1,reg_nlo,reg_nnlo,reg_n3lo);
    
    // multiplying by universal prefactor (equal to effective theory born divided by a_s^2)
    
    reg=reg * _input->_prefactor;
    
    // multiply by C^2
    double reg_ew_used;

        reg_ew_used=0.0;
    reg = reg * _input->_wc.c() * _input->_wc.c();
    
    
    if (abs(_input->_log_mur_over_muf_sq) > 1e-15)
        // mu_r evolution
    {
        
        reg = SpecialMurEvolution(reg,-_input->_log_mur_over_muf_sq);

    }
    
    // multiply each term with a_s/pi to the proper power
    
    reg.MultiplyAs(_input->_as_over_pi);
    //truncate to the requested power
    
    reg.Truncate(_input->_int_qcd_perturbative_order);
    
    double reg_nlo_used;

        reg_nlo_used = _K*reg.term_of_order(3).val();
    
    double finalReg = reg_nlo_used
    + _K*reg.term_of_order(4).val()
    + _K*reg.term_of_order(5).val()
    + reg_ew_used;
    
    
    vector<double> res;
    vector<double> lumiAtZ =  lumi_->luminosity_vector(x1,tau_/x1/z,_muF);
    
    for (int i=0;i<_number_of_components;i++){
        res.push_back(  measure * lumiAtZ[i] * finalReg);
    }


    return res;
}

HiggsEFTQQfast::HiggsEFTQQfast(){
    set_dimensions(2);
    _channel = "qq";
    _name="fast";
    _k_computed=false;
    set_rel_accuracy_multiplier(100.0);

}

void HiggsEFTQQfast::SetK(){
    _K=_input->_top_only_LO_coefficient;
    _k_computed = true;

}


vector<double> HiggsEFTQQfast::evaluateIntegral(const double* xx){
    const double x1= xx[0];
    const double z = xx[1];
    double measure = 1./x1;
    
    
    if (not(_k_computed)) SetK();

    
    // regular terms
    
    
    
    //nnlo
    
    // L^0
    const double reg_nnlo =HEFT::qq_n2lo_lz0_L0(z)
    +HEFT::qq_n2lo_lz1_L0(z)
    +HEFT::qq_n2lo_lz2_L0(z)
    +
    (
     HEFT::qq_n2lo_lz0_L1(z)
     +HEFT::qq_n2lo_lz1_L1(z)
     ) * _Lf
    +
    (
     HEFT::qq_n2lo_lz0_L2(z)
     ) * pow(_Lf,2.);
    
    const double reg_n3lo = evaluator.n3lo_reg_complete_qq(z,0)
    -evaluator.n3lo_reg_complete_qq(z,1) * _Lf
    +evaluator.n3lo_reg_complete_qq(z,2) * pow(_Lf,2.)
    -evaluator.n3lo_reg_complete_qq(z,3) * pow(_Lf,3.);

    // qq starts at a_s^2
    AsSeries reg(2,reg_nnlo,reg_n3lo);
    
    // multiplying by universal prefactor (equal to effective theory born divided by a_s^2)
    
    reg=reg * _input->_prefactor;
    
    double reg_ew_used=0.0;
    // multiply by C^2
    reg = reg * _input->_wc.c() * _input->_wc.c();
    
    
    if (abs(_input->_log_mur_over_muf_sq) > 1e-15)
        // mu_r evolution
    {
        
        reg = SpecialMurEvolution(reg,-_input->_log_mur_over_muf_sq);

    }
    
    // multiply each term with a_s/pi to the proper power
    
    reg.MultiplyAs(_input->_as_over_pi);
    //truncate to the requested power
    
    reg.Truncate(_input->_int_qcd_perturbative_order);
    
    // multiply by eft rescaling coefficient
    reg=reg* _K;
    
    double finalReg = reg.AddUp().val()+reg_ew_used;
    
    vector<double> res;
    vector<double> lumiAtZ =  lumi_->luminosity_vector(x1,tau_/x1/z,_muF);
    
    for (int i=0;i<_number_of_components;i++){
        res.push_back(  measure * lumiAtZ[i] * finalReg);
    }

    return res;
}

HiggsEFTQ1Q2fast::HiggsEFTQ1Q2fast(){
    set_dimensions(2);
    _channel = "q1q2";
    _name="fast";
    _k_computed=false;
    set_rel_accuracy_multiplier(100.0);

}

void HiggsEFTQ1Q2fast::SetK(){
    _K=_input->_top_only_LO_coefficient;
    _k_computed = true;

}

vector<double> HiggsEFTQ1Q2fast::evaluateIntegral(const double* xx){
    const double x1= xx[0];
    const double z = xx[1];
    double measure = 1./x1;
    
    if (not(_k_computed)) SetK();
    
    
    // regular terms
    
    
    
    //nnlo
    
    // L^0
    const double reg_nnlo =HEFT::q1q2_n2lo_lz0_L0(z)
    +HEFT::q1q2_n2lo_lz1_L0(z)
    +HEFT::q1q2_n2lo_lz2_L0(z)
    +
    (
     HEFT::q1q2_n2lo_lz0_L1(z)
     +HEFT::q1q2_n2lo_lz1_L1(z)
     ) * _Lf
    +
    (
     HEFT::q1q2_n2lo_lz0_L2(z)
     ) * pow(_Lf,2.);
    
    const double reg_n3lo = evaluator.n3lo_reg_complete_q1q2(z,0)
    -evaluator.n3lo_reg_complete_q1q2(z,1) * _Lf
    +evaluator.n3lo_reg_complete_q1q2(z,2) * pow(_Lf,2.)
    -evaluator.n3lo_reg_complete_q1q2(z,3) * pow(_Lf,3.);

    // q1q2 starts at a_s^2
    AsSeries reg(2,reg_nnlo,reg_n3lo);
    
    // multiplying by universal prefactor (equal to effective theory born divided by a_s^2)
    
    reg=reg * _input->_prefactor;
    
    // multiply by C^2
    double reg_ew_used=0.0;

    reg = reg * _input->_wc.c() * _input->_wc.c();
    
    
    if (abs(_input->_log_mur_over_muf_sq) > 1e-15)
        // mu_r evolution
    {
        
        reg = SpecialMurEvolution(reg,-_input->_log_mur_over_muf_sq);

    }
    
    // multiply each term with a_s/pi to the proper power
    
    reg.MultiplyAs(_input->_as_over_pi);
    //truncate to the requested power
    
    reg.Truncate(_input->_int_qcd_perturbative_order);
    
    // multiply by eft rescaling coefficient
    reg=reg* _K;
    
    double finalReg = reg.AddUp().val()+ reg_ew_used;
    
    vector<double> res;
    vector<double> lumiAtZ =  lumi_->luminosity_vector(x1,tau_/x1/z,_muF);
    
    for (int i=0;i<_number_of_components;i++){
        res.push_back(  measure * lumiAtZ[i] * finalReg);
    }
    
    return res;}



