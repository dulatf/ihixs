
#include "src/higgs/ggf_manager/ggf_manager.h"

GgfManager::GgfManager(const UserInterface& UI){
    
    _eft_input = NULL;
    mydata.Add("mh",UI.giveDouble("m_higgs"), "silent");
    mydata.Add("Etot",UI.giveInt("Etot"), "silent");
    mydata.Add("PDF set",UI.giveString("pdf_set"), "silent");
    mydata.Add("PDF member",UI.giveInt("pdf_member"), "silent");
    mydata.Add("mur",UI.giveDouble("mur"));
    mydata.Add("muf",UI.giveDouble("muf"));
    
    single(UI);

    cout << "-------------------------------------------------" << endl;
    cout << output.str() << endl;
    cout << mydata.RestrictedOutput() << endl;
    cout << "-------------------------------------------------" << endl;
    cout << " Higgs_XS = " << _Higgs_XS.val();
    if (mydata.Exists("Total Uncertainty +")) {
    cout << setprecision(4)
        << " +" << mydata.GiveVal("Total Uncertainty +").val()
        << setprecision(2)
        << "(" << mydata.GiveVal("Total Uncertainty + %").val() << "%)"
        << setprecision(4)
        << " " << mydata.GiveVal("Total Uncertainty -").val()
        << setprecision(2)
        << "(" << mydata.GiveVal("Total Uncertainty - %").val() << "%)";
    }
    cout << endl;
    cout << messages.str() << endl;
    
    stringstream output_for_file;
    output_for_file << output.str();
    output_for_file << mydata.FullOutput() << endl;
    output_for_file << messages.str() << endl;
    
    
    write_output(output_for_file.str(),
                 UI.giveString("output_filename"),
                 UI.giveString("input_filename"),
                 mydata.str(),
                 UI  );
}

void GgfManager::single(UserInterface UI) {
    
    
    
    _perform_pdf_error = UI.giveBool("with_pdf_error");
    UI.SetOption("with_pdf_error","false");
    
    if (UI.giveBool("with_eft")) perform_eft(UI);
    _qcd_perturbative_order_int = determine_int_qcd_perturbative_order( UI.giveString("qcd_perturbative_order") );
    
    _we_do_mt_expansion = false;
    if ( UI.giveBool("with_mt_expansion")) {
        if (_qcd_perturbative_order_int > 3)
            _we_do_mt_expansion = true;
        else
            messages << "note: you've requested <with_mt_expansion> but qcd_perturbative_order was set to "
            << UI.giveString("qcd_perturbative_order") << ". Hence 1/mt terms (O(NNLO)) were not computed."
            << endl;
    }
    
    
    
    if (UI.giveBool("with_exact_qcd_corrections")) perform_qcd_exact(UI);
    if (_we_do_mt_expansion) perform_mt_expansion(UI);
    if (UI.giveBool("with_resummation") and UI.giveBool("with_eft")) perform_catani_resummation(UI);
    if (UI.giveBool("with_scet")) perform_scet(UI);

    
    //: computing total
    if (UI.giveBool("with_eft")) compute_total_xs(UI);
    if (UI.giveBool("with_eft")) compute_uncertainty(UI);

    
    
    
    
   
}

void GgfManager::compute_total_xs(const UserInterface& UI) {
    //: plain fixed order
    if (_qcd_perturbative_order_int == 2)
        _Higgs_XS = mydata.GiveVal("R_LO*eftlo");
    if (_qcd_perturbative_order_int == 3)
        _Higgs_XS = mydata.GiveVal("R_LO*eftnlo");
    if (_qcd_perturbative_order_int == 4)
        _Higgs_XS = mydata.GiveVal("R_LO*eftnnlo");
    if (_qcd_perturbative_order_int == 5)
        _Higgs_XS = mydata.GiveVal("R_LO*eftn3lo");
    if (UI.giveBool("with_scet"))
        _Higgs_XS = _Higgs_XS
                +  mydata.GiveVal("R_LO") * mydata.GiveVal("delta scet total") ;
    //: fixed order +  resummation
    if (UI.giveBool("with_resummation")) {
        _Higgs_XS = mydata.GiveVal("R_LO*(eft+resum)");
    }
    if (UI.giveBool("with_exact_qcd_corrections")) {
        _Higgs_XS = _Higgs_XS + mydata.GiveVal("NLO quark mass effects");
    }
    if (UI.giveBool("with_ew_corrections")){
        _Higgs_XS = _Higgs_XS + mydata.GiveVal("ew rescaled");
    }
    if (_we_do_mt_expansion) {
        _Higgs_XS = _Higgs_XS + mydata.GiveVal("NNLO top mass effects");
    }
    
    //: adding XS
    mydata.Add("Higgs XS",_Higgs_XS);

}

void GgfManager::compute_uncertainty(const UserInterface& UI) {
    ResultPair delta_theory_per_cent;
    
    
    if (UI.giveBool("with_exact_qcd_corrections") )
        delta_theory_per_cent = delta_theory_per_cent + compute_delta_tbc(UI);
    
    if (_we_do_mt_expansion) {
        mydata.Add("delta(1/m_t)",0.01*_Higgs_XS);
        mydata.Add("delta(1/m_t) % ",ResultPair(1.0,0.0));

        delta_theory_per_cent = delta_theory_per_cent + ResultPair(1.0,0.0) ;
    }
    
    if (UI.giveBool("with_ew_corrections") )
        delta_theory_per_cent = delta_theory_per_cent + compute_delta_ew(UI);

//    if (UI.giveBool("with_truncation_error"))
//        delta_theory_per_cent = delta_theory_per_cent + compute_delta_trunc(UI);

    if (UI.giveBool("with_delta_pdf_th"))
        delta_theory_per_cent = delta_theory_per_cent + compute_delta_pdf_th(UI);
    
    
    ResultPair delta_theory_per_cent_plus = delta_theory_per_cent;
    ResultPair delta_theory_per_cent_minus = ResultPair(0.0,0.0)-delta_theory_per_cent;
    
    //: scale variation estimate
    if (UI.giveBool("with_scale_variation") ) {
        perform_scale_variation(UI);
        delta_theory_per_cent_plus = delta_theory_per_cent_plus + mydata.GiveVal("delta(scale)+(%)");
        delta_theory_per_cent_minus = delta_theory_per_cent_minus + mydata.GiveVal("delta(scale)-(%)");
    }
    
    if (_perform_pdf_error) perform_pdf_error(UI);
    
    if (UI.giveBool("with_a_s_error")) perform_as_error(UI);
    
    
    
    //: adding uncertainties (only if there is any)
    bool theory_uncertainty_computed = (UI.giveBool("with_exact_qcd_corrections")
                                        or _we_do_mt_expansion
                                        or UI.giveBool("with_ew_corrections")
                                        or UI.giveBool("with_delta_pdf_th")
                                        or UI.giveBool("with_scale_variation")) ;
    if (theory_uncertainty_computed) {
        mydata.Add("Theory Uncertainty  +",delta_theory_per_cent_plus * _Higgs_XS * 0.01);
        mydata.Add("Theory Uncertainty  -",delta_theory_per_cent_minus * _Higgs_XS * 0.01);
        mydata.Add("Theory Uncertainty % +",delta_theory_per_cent_plus);
        mydata.Add("Theory Uncertainty % -",delta_theory_per_cent_minus);
    }
    double delta_pdf_a_s_plus_percent = 0.0;
    double delta_pdf_a_s_minus_percent = 0.0;
    if (UI.giveBool("with_a_s_error")) {
        delta_pdf_a_s_plus_percent += pow(mydata.GiveVal("delta(as)+(%)").val(),2.);
        delta_pdf_a_s_minus_percent += pow(mydata.GiveVal("delta(as)-(%)").val(),2.);
    }
    if (_perform_pdf_error) {
        delta_pdf_a_s_plus_percent += pow(mydata.GiveVal("deltaPDF+(%)").val(),2.);
        delta_pdf_a_s_minus_percent += pow(mydata.GiveVal("deltaPDF-(%)").val(),2.);
    }
    delta_pdf_a_s_plus_percent = sqrt(delta_pdf_a_s_plus_percent);
    
    delta_pdf_a_s_minus_percent = sqrt(delta_pdf_a_s_minus_percent);
    
    if (UI.giveBool("with_a_s_error") or _perform_pdf_error) {
        mydata.Add("delta(PDF+a_s) +",delta_pdf_a_s_plus_percent * _Higgs_XS*0.01 );
        mydata.Add("delta(PDF+a_s) -",-delta_pdf_a_s_minus_percent * _Higgs_XS*0.01 );
        mydata.Add("delta(PDF+a_s) + %",ResultPair(delta_pdf_a_s_plus_percent,0.0));
        mydata.Add("delta(PDF+a_s) - %",ResultPair(-delta_pdf_a_s_minus_percent,0.0));
        if (theory_uncertainty_computed) {
            mydata.Add("Total Uncertainty +",delta_pdf_a_s_plus_percent * _Higgs_XS*0.01 + delta_theory_per_cent_plus * _Higgs_XS*0.01 );
            mydata.Add("Total Uncertainty -",-delta_pdf_a_s_minus_percent * _Higgs_XS*0.01 + delta_theory_per_cent_minus * _Higgs_XS*0.01);
            mydata.Add("Total Uncertainty + %",ResultPair(delta_pdf_a_s_plus_percent,0.0) + delta_theory_per_cent_plus);
            mydata.Add("Total Uncertainty - %",ResultPair(-delta_pdf_a_s_minus_percent,0.0) + delta_theory_per_cent_minus);
        }
    }
    else {
        if (theory_uncertainty_computed) {
            mydata.Add("Total Uncertainty +", delta_theory_per_cent_plus * _Higgs_XS*0.01 );
            mydata.Add("Total Uncertainty -", delta_theory_per_cent_minus * _Higgs_XS*0.01);
            mydata.Add("Total Uncertainty + %", delta_theory_per_cent_plus);
            mydata.Add("Total Uncertainty - %", delta_theory_per_cent_minus);
        }
    }
        
    
}

ResultPair GgfManager::compute_delta_tbc(const UserInterface& UI) {
    if (_we_do_mt_expansion) {
        ResultPair delta_tbc = mydata.GiveVal("delta tbc ratio")
            * (  mydata.GiveVal("R_LO*eftnnlo")
               - mydata.GiveVal("R_LO*eftnlo")
               + mydata.GiveVal("NNLO top mass effects")
               );
        // the magic_factor accounts for scheme dependence
        const double magic_factor = 1.3;
        mydata.Add("delta_tbc",delta_tbc * magic_factor);
        mydata.Add("delta_tbc %",
                   delta_tbc/_Higgs_XS * 100.0 * magic_factor);
        return  delta_tbc/_Higgs_XS * 100.0 * magic_factor;
    }
    else {
        messages << "note: delta_tbc was requested (an estimate of the uncertainty due to the lack of light quark mass effects at NNLO) but either the with_mt_expansion is set to false or the qcd perturbative order is less than NNLO. In either case we cannot estimate delta_tbc, so it is set to zero" << endl;
        return ResultPair(0.0,0.0);
    }

}

ResultPair GgfManager::compute_delta_ew(const UserInterface& UI) {
    ResultPair delta_EW_per_cent(1.0,0.0);
    ResultPair delta_EW_absolute = delta_EW_per_cent * _Higgs_XS * 0.01;

    mydata.Add("delta EW",delta_EW_absolute);
    mydata.Add("delta EW %",delta_EW_per_cent);
    return delta_EW_per_cent;


}


ResultPair GgfManager::compute_delta_pdf_th(const UserInterface& UI) {
    if (UI.giveString("qcd_perturbative_order")  == "N3LO" ) {
        perform_pdf_th_error(UI);
        ResultPair delta_pdf_th =  (
                                    mydata.GiveVal("R_LO*eftnnlo (with NLO PDF)")
                                    - mydata.GiveVal("R_LO*eftnnlo")
                                    )
        / mydata.GiveVal("R_LO*eftnnlo")
        *  0.5 * 100.;
        
        delta_pdf_th.TakeAbsoluteValue();
        
        mydata.Add("delta PDF-TH %",delta_pdf_th);
        
        return  delta_pdf_th;
    }
    else {
        messages << "note: with_delta_pdf_theory was requested "
        << "but qcd_perturbative_order = " << UI.giveString("qcd_perturbative_order")
        << ", so PDF-THEORY error was not computed." << endl;
        return ResultPair(0.0,0.0);
        }
    
}

void GgfManager::perform_eft(const UserInterface& UI){
    //: EFT
    InclusiveHiggsEFT eft(UI);
    //: printing input information
    //output<<eft.InputInformation();
    //cout << "as(muf) = " << eft.AsOverPi()*consts::Pi << endl;
    cout << "Computing " <<UI.giveString("qcd_perturbative_order")<< " cross section in the Effective Field Theory approach" <<endl;
    eft.Evaluate();
    
    if (UI.giveBool("with_eft_channel_info")) {cout << eft << endl;}
    
    mydata.Add("as_at_mz",eft.GetInput()->_as_at_mz, "silent");
    mydata.Add("as_at_mur",eft.GetInput()->_as_over_pi * consts::Pi, "silent");
    mydata.Add("mt_used",eft.GetInput()->_model.top.m(), "silent");
    mydata.Add("mb_used",eft.GetInput()->_model.bottom.m(), "silent");
    mydata.Add("mc_used",eft.GetInput()->_model.charm.m(), "silent");

    AsSeries eft_res =eft.QCDCorrections();
    AsSeries eft_ew =eft.EwCorrections();
    
    
    double K =eft.RescalingCoefficient();
    
    
    
    ResultPair eft_lo = eft_res.term_of_order(2);
    ResultPair eft_n1lo = eft_res.term_of_order(2)+eft_res.term_of_order(3);
    ResultPair eft_n2lo = eft_res.term_of_order(2)+eft_res.term_of_order(3)+eft_res.term_of_order(4);
    ResultPair eft_n3lo = eft_res.term_of_order(2)+eft_res.term_of_order(3)+eft_res.term_of_order(4)+eft_res.term_of_order(5);
    ResultPair eft_gg = eft.GiveChannel("gg")->QCDResult().AddUp();
    ResultPair eft_qg = eft.GiveChannel("qg")->QCDResult().AddUp();
    ResultPair eft_qqbar = eft.GiveChannel("qqbar")->QCDResult().AddUp();
    ResultPair eft_qq = eft.GiveChannel("qq")->QCDResult().AddUp();
    ResultPair eft_q1q2 = eft.GiveChannel("q1q2")->QCDResult().AddUp();
    
    
    int qcd_perturbative_order_int = determine_int_qcd_perturbative_order(UI.giveString("qcd_perturbative_order"));
    
    mydata.Add("eftlo",eft_lo, "silent");
    if (qcd_perturbative_order_int > 2) mydata.Add("eftnlo",eft_n1lo, "silent");
    if (qcd_perturbative_order_int > 3) mydata.Add("eftnnlo",eft_n2lo, "silent");
    if (qcd_perturbative_order_int > 4) mydata.Add("eftn3lo",eft_n3lo, "silent");
    mydata.Add("R_LO",K);
    mydata.Add("R_LO*eftlo", K * eft_lo);
    if (qcd_perturbative_order_int > 2) mydata.Add("R_LO*eftnlo", K * eft_n1lo);
    if (qcd_perturbative_order_int > 3) mydata.Add("R_LO*eftnnlo", K * eft_n2lo);
    if (qcd_perturbative_order_int > 4) mydata.Add("R_LO*eftn3lo", K * eft_n3lo);
    if (qcd_perturbative_order_int > 2) mydata.Add("ggnlo/eftnlo",(eft.GiveChannel("gg")->QCDResult().term_of_order(2)+
                                                                   eft.GiveChannel("gg")->QCDResult().term_of_order(3))/eft_n1lo, "silent");
    if (qcd_perturbative_order_int > 2) mydata.Add("qgnlo/eftnlo",(eft.GiveChannel("qg")->QCDResult().term_of_order(2)+
                                                                   eft.GiveChannel("qg")->QCDResult().term_of_order(3))/eft_n1lo, "silent");
    
    
    if (qcd_perturbative_order_int > 3) mydata.Add("ggnnlo/eftn2lo",(eft.GiveChannel("gg")->QCDResult().term_of_order(2)+
                                                                     eft.GiveChannel("gg")->QCDResult().term_of_order(3)+
                                                                     eft.GiveChannel("gg")->QCDResult().term_of_order(4))/eft_n2lo, "silent");
    if (qcd_perturbative_order_int > 3) mydata.Add("qgnnlo/eftn2lo",(eft.GiveChannel("qg")->QCDResult().term_of_order(2)+
                                                                     eft.GiveChannel("qg")->QCDResult().term_of_order(3)+
                                                                     eft.GiveChannel("qg")->QCDResult().term_of_order(4))/eft_n2lo, "silent");
    if (qcd_perturbative_order_int > 4) mydata.Add("ggn3lo/eftn3lo",eft_gg/eft_n3lo, "silent");
    if (qcd_perturbative_order_int > 4) mydata.Add("qgn3lo/eftn3lo",eft_qg/eft_n3lo, "silent");
    
    // necessary for mt_expansion
    if (UI.giveBool("with_mt_expansion") and qcd_perturbative_order_int > 3) {
        mydata.Add("gg at order 4",
                   K * eft.GiveChannel("gg")->QCDResult().term_of_order(4),
                   "internal"
                   );
        mydata.Add("qg at order 4",
                   K * eft.GiveChannel("qg")->QCDResult().term_of_order(4),
                   "internal"
                   );
    }
    mydata.Add("R_LO*gg channel", K * eft_gg,"silent");
    mydata.Add("R_LO*qg channel",K * eft_qg, "silent");
    mydata.Add("R_LO*qqbar channel",K * eft_qqbar, "silent");
    mydata.Add("R_LO*qq channel",K * eft_qq, "silent");
    mydata.Add("R_LO*q1q2 channel",K * eft_q1q2, "silent");
    
    //: ew corrections
    if (UI.giveBool("with_ew_corrections")){
        mydata.Add("ew rescaled as^2",eft_ew.term_of_order(2)*K, "silent");
        if (qcd_perturbative_order_int > 2) mydata.Add("ew rescaled as^3",eft_ew.term_of_order(3)*K, "silent");
        if (qcd_perturbative_order_int > 3) mydata.Add("ew rescaled as^4",eft_ew.term_of_order(4)*K, "silent");
        if (qcd_perturbative_order_int > 4) mydata.Add("ew rescaled as^5",eft_ew.term_of_order(5)*K, "silent");
        
        if (qcd_perturbative_order_int > 2) mydata.Add("mixed EW-QCD",
                                                       K * ( eft_ew.AddUp()-eft_ew.term_of_order(2) ), "silent"
                                                       );
        mydata.Add("ew rescaled",eft_ew.AddUp()*K);
        
        if (qcd_perturbative_order_int > 2) {
            ResultPair d_s_gg_NLO_reg=eft.GiveChannel("gg")->TermPtr("gg_reg_nlo")->QCDResult().term_of_order(3);
            ResultPair d_s_qg_NLO_reg=eft.GiveChannel("qg")->TermPtr("qg_reg_nlo")->QCDResult().term_of_order(3);
            ResultPair d_s_qqbar_NLO_reg=eft.GiveChannel("qqbar")->TermPtr("qqb_reg_nlo")->QCDResult().term_of_order(3);
            ResultPair d_s_NLO = eft_res.term_of_order(3);
            
            ResultPair hard_ratio =  (
                                      d_s_gg_NLO_reg +
                                      d_s_qg_NLO_reg +
                                      d_s_qqbar_NLO_reg
                                      )
            / d_s_NLO;
            if ( hard_ratio.val() < 0 ) hard_ratio=ResultPair(-hard_ratio.val(),hard_ratio.err());
            mydata.Add("hard ratio from eft", hard_ratio, "silent");
        }
    }
    //: factorized form
    //WC.Truncate(5);
    //ResultPair WilsonCoeff = WC.AddUp();
    AsSeries WC = eft.WC();
    //cout << "WC= " << WC << endl;
    WC.Truncate(qcd_perturbative_order_int-2);
    
    
    AsSeries Eta = eft.Eta();
    //cout << "n= " << Eta << endl;
    Eta.Truncate(qcd_perturbative_order_int);
    
    //AsSeries eft_recalc = WC*WC*Eta;
    //eft_recalc.Truncate(5);
    
    
    AsSeries eft_factorized = WC*WC*Eta;
    ResultPair eft_sigma_fac = eft_factorized.AddUp();
    mydata.Add("WC",WC.AddUp(), "silent");
    mydata.Add("WC^2",WC.AddUp()*WC.AddUp(), "silent");
    AsSeries WCsquared_trunc = WC*WC;
    WCsquared_trunc.Truncate(qcd_perturbative_order_int-2);
    //cout << "|C|^2 = " << WCsquared_trunc << endl;
    mydata.Add("WC^2_trunc",WCsquared_trunc.AddUp(), "silent");
    mydata.Add("n",Eta.AddUp(), "silent");
    //mydata.Add("Cbern",pow(10,4.0) * pow( eft.AsOverPi() / 3 , 2.0) * WCsquared_trunc.AddUp() , "silent"    );
    mydata.Add("sigma factorized",eft_sigma_fac, "silent");
    //mydata.Add("sigma recalc",eft_recalc.AddUp());
    
    // printing out detailed view per integral
    if (UI.giveString("verbose") == "medium")
        output << eft << endl;
    
    //IMPORTANT: setting _eft_input which is necessary for resummation
    _eft_input = eft.GetInput();
    
    
    
}

vector<double> GgfManager::compute_scale_variation_boundary(const string order) {
    vector<double> boundary;
    if (order == "N3LO") boundary.push_back(0.4);
    else boundary.push_back(0.25);
    
    boundary.push_back(1.0);
    return boundary;
}

string GgfManager::make_string_for_scale_variation_boundary(const double& val) {
    std::ostringstream s; s << val << "*mh";
    return s.str();
}


void GgfManager::perform_scale_variation(const UserInterface& UI) {
    
    double mur_user = UI.giveDouble("mur");
    double muf_user = UI.giveDouble("muf");
    double mh_central = UI.giveDouble("m_higgs");
    
    
    if ( mur_user == mh_central / 2.0 and muf_user == mh_central / 2.0) {
        vector<string> name_of_maximal_eft_order_rescaled = {"R_LO*eftlo",
            "R_LO*eftnlo",
            "R_LO*eftnnlo",
            "R_LO*eftn3lo"};
        
        int int_pert_order = determine_int_qcd_perturbative_order(UI.giveString("qcd_perturbative_order"))-2;
        ResultPair central = mydata.GiveVal(name_of_maximal_eft_order_rescaled[int_pert_order]);
        
        
        vector<double> boundary = compute_scale_variation_boundary(UI.giveString("qcd_perturbative_order"));
        UserInterface UI_low_scale = UI;
        UI_low_scale.SetOption("mur",  boundary[0] * mh_central);
        UI_low_scale.SetOption("muf",  boundary[0] * mh_central);
        InclusiveHiggsEFT* eft_low = new InclusiveHiggsEFT(UI_low_scale);
        cout << "Scale variation: Computing " << UI.giveString("qcd_perturbative_order")
        << " cross section at mu = "
        << boundary[0] << " * mh" << endl;
        eft_low->Evaluate();
        ResultPair low_scale = eft_low->R_LO_x_EFT_N3LO();
        
        mydata.Add("rEFT(low)", low_scale, "silent");
        
        
        double low_end = low_scale.val() - central.val();
        double low_end_percent = low_end / central.val() * 100.0;
        
        
        
        //computing eft at mu =  mh
        UserInterface UI_high_scale = UI;
        UI_high_scale.SetOption("mur" ,  mh_central * boundary[1]);
        UI_high_scale.SetOption("muf" ,  mh_central * boundary[1]);
        InclusiveHiggsEFT* eft_high = new InclusiveHiggsEFT(UI_high_scale);
        cout << "Scale variation: Computing " << UI.giveString("qcd_perturbative_order")
        << " cross section at mu = " << boundary[1]
        << " mh" << endl;
        eft_high->Evaluate();
        ResultPair high_scale = eft_high->R_LO_x_EFT_N3LO();
        mydata.Add("rEFT(high)", high_scale, "silent");
        double high_end = high_scale.val() - central.val();
        double high_end_percent = high_end / central.val() * 100.0;
        
        mydata.Add("delta(scale)+",low_end);
        mydata.Add("delta(scale)-",high_end);
        mydata.Add("delta(scale)+(%)",low_end_percent);
        mydata.Add("delta(scale)-(%)",high_end_percent);
        
        compute_scale_variation_for_unrescaled_eft(eft_low,eft_high,int_pert_order);
        //cout << "with_lower_ord_scale_var = " << UI.giveString("with_lower_ord_scale_var") << endl;
        if (UI.giveBool("with_lower_ord_scale_var")) {
            InclusiveHiggsEFT* eft_low_special;
            if (int_pert_order == 3) {
                UserInterface curUI = UI;
                curUI.SetOption("mur" ,  mh_central * 0.25);
                curUI.SetOption("muf" ,  mh_central * 0.25);
                curUI.SetOption("qcd_perturbative_order" ,  "NNLO");
                cout << "Lower order scale variation: Computing " << curUI.giveString("qcd_perturbative_order")
                << " cross section at mu = " << 0.25 << " mh" << endl;
                eft_low_special = new InclusiveHiggsEFT(curUI);
                eft_low_special->Evaluate();
            }
            else {
                eft_low_special = eft_low;
            }
            for (int i=0; i < int_pert_order; i++) {
                ResultPair cur_high_scale = eft_high->R_LO_x_EFT(i+2);
                ResultPair cur_low_scale = eft_low_special->R_LO_x_EFT(i+2);
                ResultPair cur_central = mydata.GiveVal(name_of_maximal_eft_order_rescaled[i]);
                double cur_high_end = cur_high_scale.val() - cur_central.val();
                double cur_high_end_percent = cur_high_end / cur_central.val() * 100.0;
                double cur_low_end = cur_low_scale.val() - cur_central.val();
                double cur_low_end_percent = cur_low_end / cur_central.val() * 100.0;
                vector<string> ordername = {"LO","NLO","NNLO"};
                mydata.Add(string("delta(scale)+at")+ordername[i],cur_low_end,"silent");
                mydata.Add(string("delta(scale)-at")+ordername[i],cur_high_end,"silent");
                mydata.Add(string("delta(scale)+(%)at")+ordername[i],cur_low_end_percent,"silent");
                mydata.Add(string("delta(scale)-(%)at")+ordername[i],cur_high_end_percent,"silent");
            }
        }
        
    }
    else {
        messages << "The requested scale variation estimate was not computed";
        if (mur_user != mh_central / 2.)
            messages << " because mur = " << mur_user
            << " != "
            << mh_central/2. << " = m_h/2 " << endl;
        else
            messages << " because muf = " << muf_user
            << " != "
            << mh_central/2. << " = m_h/2 " << endl;
        
    }
}

void GgfManager::compute_scale_variation_for_unrescaled_eft(InclusiveHiggsEFT* eft_low, InclusiveHiggsEFT* eft_high,int int_pert_order){
    // computing scale uncertainty of pure eft
    // retrieving central
    vector<string> name_of_maximal_eft_order = {"eftlo",
        "eftnlo",
        "eftnnlo",
        "eftn3lo"};
    
    ResultPair central_pure_eft = mydata.GiveVal(name_of_maximal_eft_order[int_pert_order]);
    ResultPair low_scale_pure_eft = eft_low->EFT_N3LO();
    ResultPair high_scale_pure_eft = eft_high->EFT_N3LO();
    double high_end_pure_eft = high_scale_pure_eft.val() - central_pure_eft.val();
    double high_end_pure_eft_percent = high_end_pure_eft / central_pure_eft.val() * 100.0;
    double low_end_pure_eft = low_scale_pure_eft.val() - central_pure_eft.val();
    double low_end_pure_eft_percent = low_end_pure_eft / central_pure_eft.val() * 100.0;

    mydata.Add("delta(scale)+ pure eft",low_end_pure_eft,"silent");
    mydata.Add("delta(scale)- pure eft",high_end_pure_eft,"silent");
    mydata.Add("delta(scale)+(%) pure eft",low_end_pure_eft_percent,"silent");
    mydata.Add("delta(scale)-(%) pure eft",high_end_pure_eft_percent,"silent");
}

void GgfManager::perform_as_error(const UserInterface& UI) {
    
    bool as_error_available_for_this_pdf = false;
    UserInterface UIlocalplus = UI;
    UserInterface UIlocalminus = UI;

    if (UI.giveString("pdf_set")=="PDF4LHC15_nnlo_100_pdfas") {
        UIlocalplus.SetOption("pdf_member",102);
        UIlocalminus.SetOption("pdf_member",101);
        as_error_available_for_this_pdf = true;
    }
    if (UI.giveString("pdf_set")=="PDF4LHC15_nnlo_100") {
        UIlocalplus.SetOption("pdf_set","PDF4LHC15_nnlo_100_pdfas");
        UIlocalplus.SetOption("pdf_member",102);
        UIlocalminus.SetOption("pdf_set","PDF4LHC15_nnlo_100_pdfas");
        UIlocalminus.SetOption("pdf_member",101);
        as_error_available_for_this_pdf = true;
    }
    
    if (as_error_available_for_this_pdf) {
        InclusiveHiggsEFT* eft_asplus = new InclusiveHiggsEFT(UIlocalplus);
        cout << "Computing eft cross section with a_s + delta(a_s) " << endl;
        eft_asplus->Evaluate();
        ResultPair xs_plus = eft_asplus->R_LO_x_EFT_N3LO();
        mydata.Add("rEFT(as+)", xs_plus,"silent");

        
        InclusiveHiggsEFT* eft_asminus = new InclusiveHiggsEFT(UIlocalminus);
        cout << "Computing eft cross section with a_s - delta(a_s) " << endl;
        eft_asminus->Evaluate();
        ResultPair xs_minus = eft_asminus->R_LO_x_EFT_N3LO();
        mydata.Add("rEFT(as-)", xs_minus,"silent");
        
        
        vector<string> name_of_maximal_eft_order_rescaled = {"R_LO*eftlo",
            "R_LO*eftnlo",
            "R_LO*eftnnlo",
            "R_LO*eftn3lo"};
        
        int int_pert_order = determine_int_qcd_perturbative_order(UI.giveString("qcd_perturbative_order"))-2;
        ResultPair central = mydata.GiveVal(name_of_maximal_eft_order_rescaled[int_pert_order]);
        
        double deltaplus = xs_plus.val() - central.val();
        double deltaminus = -central.val() + xs_minus.val();
        mydata.Add("delta(as)+",deltaplus);
        mydata.Add("delta(as)-",deltaminus);
        mydata.Add("delta(as)+(%)",deltaplus/ central.val() * 100.0);
        mydata.Add("delta(as)-(%)",deltaminus/ central.val() * 100.0);
        
    }
    else {
        cout << "a_s uncertainty unavailable for " << UI.giveString("pdf_set")
            << endl;
        messages << "The requested a_s uncertainty estimate was not computed: ";
        messages << "this feature is not available for the pdf set "
                << UI.giveString("pdf_set") << endl;
    }
}

void GgfManager::perform_mt_expansion(const UserInterface& UI) {
    //: MtExpansion
    cout << "Computing NNLO 1/mt terms " <<endl;
    InclusiveHiggsMtExpansion* mt_expansion_run = new InclusiveHiggsMtExpansion(UI);
    mt_expansion_run->Evaluate();
    //AsSeries mt_expansion=mt_expansion_run->QCDCorrections();
    
    ResultPair DeltaMtExpansiongg =
    mt_expansion_run->delta_mt_gluon_gluon(
                                           mydata.GiveVal("gg at order 4")
                                           );
    mydata.Add("NNLO mt exp gg",DeltaMtExpansiongg, "silent");
    //output<<*mt_expansion_run<<endl;
    
    ResultPair DeltaMtExpansionqg = mt_expansion_run->delta_mt_quark_gluon(
                                                                           mydata.GiveVal("qg at order 4")
                                                                           );
    mydata.Add("NNLO mt exp qg",DeltaMtExpansionqg, "silent");
    ResultPair DeltaMtExpansionTotal = DeltaMtExpansiongg + DeltaMtExpansionqg;
    mydata.Add("NNLO top mass effects",DeltaMtExpansionTotal);
    
}

void GgfManager::perform_qcd_exact(const UserInterface& UI) {
    
    UserInterface UInlo = UI;
    //: if you ever want to compute deltaQCD with NLO pdfs uncomment below
    //if (UInlo.pdf_set_for_nlo=="none")
    //    UInlo.pdf_set_for_nlo = UInlo.pdf_set;
    //UInlo.pdf_set = UInlo.pdf_set_for_nlo;
    // Question: what perturbative order we run a_s through when computing Delta_QCD ?
    // Answer presently:  default N3LO, defined by qcd_order_evol
    //
    if (UInlo.giveString("qcd_perturbative_order")!="LO")
        UInlo.SetOption("qcd_perturbative_order",  "NLO");
    
    //: t+b+c
    perform_exact_nlo_with_selected_quarks(UInlo,UI.giveDouble("y_top"), UI.giveDouble("y_bot"), UI.giveDouble("y_charm"));
    
    //: t
    perform_exact_nlo_with_selected_quarks(UInlo,UI.giveDouble("y_top"), 0.0, 0.0);
    
    //: computing NLO quark mass effects and the uncertainty on that
    //: we only do this if the eft is already computed
    if (mydata.Exists("R_LO*eftlo")) {
        if (UInlo.giveString("qcd_perturbative_order") == "LO") {
            ResultPair deltaQCD=mydata.GiveVal("exact LO t+b+c") - mydata.GiveVal("R_LO*eftlo");
            mydata.Add("NLO quark mass effects",deltaQCD);
            mydata.Add("NLO quark mass effects / eft %",deltaQCD.val()/(mydata.GiveVal("R_LO*eftlo").val())*100.0, "silent");
        
        }
    
        if (UInlo.giveString("qcd_perturbative_order") == "NLO") {
            ResultPair deltaQCD=mydata.GiveVal("exact NLO t+b+c") - mydata.GiveVal("R_LO*eftnlo");
            ResultPair delta_sigma_tbc_NLO = mydata.GiveVal("exact NLO t+b+c") - mydata.GiveVal("exact LO t+b+c");
            mydata.Add("NLO quark mass effects",deltaQCD);
            mydata.Add("NLO quark mass effects / eft %",deltaQCD.val()/(mydata.GiveVal("R_LO*eftnlo").val())*100.0, "silent");
            // computing necessary ratio for delta_tbc
            // NLO contribution to sigma_t
            ResultPair delta_sigma_t_NLO=mydata.GiveVal("exact NLO t") - mydata.GiveVal("exact LO t");
            // ratio = | ds_{t}^NLO - ds_{tbc}^NLO) / ds_{t}^NLO |
            ResultPair d_tbc_ratio = (
                                      delta_sigma_t_NLO -
                                      delta_sigma_tbc_NLO
                                      ) / delta_sigma_t_NLO;
            if (d_tbc_ratio.val()<0) d_tbc_ratio = ResultPair(-d_tbc_ratio.val(),d_tbc_ratio.err());
            
            mydata.Add("delta sigma t NLO",delta_sigma_t_NLO, "internal");
            mydata.Add("delta sigma t+b+c NLO",delta_sigma_tbc_NLO, "internal");
            mydata.Add("delta tbc ratio",d_tbc_ratio, "internal");
        
            }
        }
    
    if (UInlo.giveBool("with_indiv_mass_effects")) {
        perform_qcd_exact_study(UInlo);
    }
}



void GgfManager::perform_qcd_exact_study(const UserInterface& UI) {
    //: t+b
    perform_exact_nlo_with_selected_quarks(UI,1.0,1.0,0.0);
    //: t+c
    perform_exact_nlo_with_selected_quarks(UI,1.0,0.0,1.0);
    //: b
    perform_exact_nlo_with_selected_quarks(UI,0.0,1.0,0.0);
    //: c
    perform_exact_nlo_with_selected_quarks(UI,0.0,0.0,1.0);
    //: b+c
    perform_exact_nlo_with_selected_quarks(UI,0.0,1.0,1.0);
}

void GgfManager::perform_exact_nlo_with_selected_quarks(const UserInterface& UI, double yt, double yb, double yc) {
    UserInterface UIlocal = UI;
    UIlocal.SetOption("y_top" ,  yt ) ;
    UIlocal.SetOption("y_bot" , yb );
    UIlocal.SetOption("y_charm" , yc );
    InclusiveHiggsExact* exact_nlo = new InclusiveHiggsExact(UIlocal);
    
    vector<string> quark_names = {"t", "b", "c"};
    vector<double> y_values = {yt, yb, yc};
    stringstream run_name;
    bool started = false;
    for (int i = 0 ; i < 3 ; i++ ) {
        if (y_values[i] > 0.0 ) {
            if (started) run_name << "+";
            run_name << quark_names[i];
            started = true;
        }
    }
    
    cout << "Computing " << UI.giveString("qcd_perturbative_order") << " exact " <<run_name.str() << endl;
    exact_nlo->Evaluate();
    
    
    string verbose_level = "silent";
    if (UI.giveBool("with_indiv_mass_effects")) verbose_level = "all";
    
    stringstream lo_name;
    lo_name << "exact LO " << run_name.str();
    mydata.Add(lo_name.str(),
               exact_nlo->QCDCorrections().term_of_order(2), verbose_level);
    if (UI.giveString("qcd_perturbative_order") == "NLO") {
        stringstream nlo_name;
        nlo_name << "exact NLO " << run_name.str();
        mydata.Add(nlo_name.str(),
                   exact_nlo->QCDCorrections().AddUp(),verbose_level);
    }
    
}




void GgfManager::perform_catani_resummation(const UserInterface& UI){
    cout << "Computing resummation contributions, flavour: " << UI.giveString("resummation_type") <<endl;
    int qcd_perturbative_order_int = determine_int_qcd_perturbative_order(UI.giveString("qcd_perturbative_order"));
    
    
    const string grid_directory="./grids/";
    if (_eft_input == NULL) {
        cout << "Error: attempted to run resummation without previously having run eft." << endl;
        exit(0);
    }
    //: threshold resummation
    ThresholdResummation trs(UI,grid_directory,_eft_input);
    //unsigned int log_order = 3; //: 2 is NNLL
    //unsigned int matchingOrder = 3;//: matching to NNLO fixed order
    bool pi2resummation = false;
    ResultPair resumation_N3LL_correction(
                    trs.ResummationCorrection(UI.giveInt("resummation_log_order"),
                    UI.giveInt("resummation_matching_order"),
                    pi2resummation)
                ,1e-6);
    
    ThresholdResummation my_resummation(UI,grid_directory,_eft_input);
    ResultPair deltaResLO   = ResultPair(trs.ResummationCorrection(0,0,false),0.0);
    ResultPair deltaResN1LO = ResultPair(trs.ResummationCorrection(1,1,false),0.0);
    ResultPair deltaResN2LO = ResultPair(trs.ResummationCorrection(2,2,false),0.0);
    ResultPair deltaResN3LO = ResultPair(trs.ResummationCorrection(3,3,false),0.0);
    ResultPair eft_LOpLL =     mydata.GiveVal("eftlo")   + deltaResLO;
    mydata.Add("lopll",eft_LOpLL);
    ResultPair rescaled_total_resummed = mydata.GiveVal("R_LO") * eft_LOpLL;
    if (qcd_perturbative_order_int > 2) {
        ResultPair eft_NLOpNLL =   mydata.GiveVal("eftnlo") + deltaResN1LO;
        mydata.Add("nlopnll",eft_NLOpNLL);
        rescaled_total_resummed = mydata.GiveVal("R_LO") * eft_NLOpNLL;
    }
    if (qcd_perturbative_order_int > 3) {
        ResultPair eft_NNLOpNNLL = mydata.GiveVal("eftnnlo") + deltaResN2LO;
        mydata.Add("nnlopnnll",eft_NNLOpNNLL);
        rescaled_total_resummed = mydata.GiveVal("R_LO") * eft_NNLOpNNLL;
        
    }
    if (qcd_perturbative_order_int > 4) {
        ResultPair eft_N3LOpN3LL = mydata.GiveVal("eftn3lo") + deltaResN3LO;
        mydata.Add("n3lopn3ll",eft_N3LOpN3LL);
        rescaled_total_resummed = mydata.GiveVal("R_LO") * eft_N3LOpN3LL;
        
    }
    mydata.Add("R_LO*(eft+resum)", rescaled_total_resummed,"silent");

    
  

}

void GgfManager::perform_pdf_error(UserInterface UI){

    UI.SetOption("with_pdf_error","true");
    UI.SetOption("with_exact_qcd_corrections","false");
    UI.SetOption("with_ew_corrections","false");
    if (UI.giveString("pdf_set")=="PDF4LHC15_nnlo_100_pdfas") {
        UI.SetOption("pdf_set","PDF4LHC15_nnlo_100");
    }


    const LHAPDF::PDFSet mypdfset(UI.giveString("pdf_set"));
    cout << "Computing PDF error with "
        << mypdfset.name() << " which has " << mypdfset.size()
        << " members" << endl
        << "Description : " << mypdfset.description() <<endl
        << "Data version : " << mypdfset.dataversion() << endl
        << "Uncertainty type : " << mypdfset.errorType()
        << endl;
    
    AsSeries res;
    InclusiveHiggsEFTFast* fast_run = new InclusiveHiggsEFTFast(UI);
    fast_run->Evaluate();
    res = fast_run->QCDCorrections();
    vector<double> results = fast_run->ResultVector();
    vector<double> errors  = fast_run->ErrorVector();
    

    
    ResultPair central_xs(results[0],errors[0]);
    
    
    const LHAPDF::PDFUncertainty delta_pdf = mypdfset.uncertainty(results, 68); // -1 => same C.L. as set, 68: 68%CL
    mydata.Add("R_LO*eftn3lo_central",central_xs,"silent");
    mydata.Add("deltaPDF+",delta_pdf.errplus);
    mydata.Add("deltaPDF-",delta_pdf.errminus);
    mydata.Add("deltaPDFsymm",delta_pdf.errsymm);
    mydata.Add("deltaPDF+(%)",delta_pdf.errplus/delta_pdf.central*100.);
    mydata.Add("deltaPDF-(%)",delta_pdf.errminus/delta_pdf.central*100.);
    
    if (UI.giveString("verbose") == "medium") {
        output<<fast_run->InputInformation();
    
        for (int i=0;i<results.size();i++){
            output<<i<<" : "<<results[i]<<" ["<<errors[i]<<"]"<<endl;
        }
    }
}



void GgfManager::perform_pdf_th_error(const UserInterface& UI) {
    
    
    if ( UI.giveString("pdf_set")==UI.giveString("pdf_set_for_nlo") ) {
        cout << "warning: pdf_set_for_nlo == pdf_set. The delta(PDF-TH) cannot be computed." << endl;
        mydata.Add("R_LO*eftnnlo (with NLO PDF)",0.0,"silent");
        return;
    }
    
    cout << "Computing estimate for PDF-TH error" <<endl;
    UserInterface UIpdfth=UI;
    UIpdfth.SetOption("pdf_set",UI.giveString("pdf_set_for_nlo"));
    UIpdfth.SetOption("qcd_perturbative_order" , "N3LO");
    UIpdfth.SetOption("with_fixed_as_at_mz", mydata.GiveVal("as_at_mz").val());
    InclusiveHiggsEFT eft(UIpdfth);
    eft.Evaluate();
    
    AsSeries eft_res =eft.QCDCorrections();
    double R_LO =eft.RescalingCoefficient();
    
    ResultPair eft_n2lo = eft_res.term_of_order(2)+eft_res.term_of_order(3)+eft_res.term_of_order(4);
    mydata.Add( "R_LO*eftnnlo (with NLO PDF)" , R_LO * eft_n2lo, "silent" );
    
    
    
    AsSeries WC = eft.WC();
    AsSeries Eta = eft.Eta();
    
    
    AsSeries factorized = WC*WC*Eta;
    

    if (UI.giveString("verbose") == "medium") output << eft << endl;
}

void GgfManager::perform_scet(const UserInterface& UI){
    
    
    cout << "Computing SCET: we use on-shell top mass = "
        << UI.giveDouble("mt_on_shell") << endl;
    const double tau = pow(UI.giveDouble("m_higgs") / UI.giveDouble("Etot") , 2.);
    const double Gf=1.16637/100000.0;
    const double MZ=91.1876;
    
    
    
    SCETTER SCET;
    
    SCET.InitiatePDF(UI.giveString("pdf_set"));
    double alphamz=SCET.GetAlpha(MZ);
    //cout << "a_s(" << MZ << ") = " << alphamz << endl;
    
    SCET.SetAlphaOrder(3,alphamz,MZ);
    cout<<setprecision(5)<<std::fixed;
    
    vector<double> deltaSCET=SCET.DeltaSCET(UI.giveDouble("m_higgs"),UI.giveDouble("mt_on_shell"),Gf,UI.giveDouble("muf"),tau);
    
    int qcd_perturbative_order_int = determine_int_qcd_perturbative_order(UI.giveString("qcd_perturbative_order"));
    mydata.Add("delta scet LO",deltaSCET[0]);
    
    if (qcd_perturbative_order_int>2)
        mydata.Add("delta scet NLO",deltaSCET[1]);
    if (qcd_perturbative_order_int>3)
        mydata.Add("delta scet NNLO",deltaSCET[2]);
    if (qcd_perturbative_order_int>4)
        mydata.Add("delta scet N3LO",deltaSCET[3]);
    mydata.Add("delta scet total",deltaSCET[qcd_perturbative_order_int-2],"silent");
}

int GgfManager::determine_int_qcd_perturbative_order(const string & order) {
    if (order=="LO")  return 2;
    if (order=="NLO") return 3;
    if (order=="NNLO") return 4;
    if (order=="N3LO") return 5;
    cout << "Error in " << __func__
    << " qcd_perturbative_order was declared wrongly : <"
    << order << ">" << endl;
    exit (EXIT_FAILURE);
    
}


string GgfManager::output_text_separator() {
    return "------------------------------------------------------------------";
}


void GgfManager::write_output(const string& output,
                              const string& filename,
                              const string& input_filename,
                              const string& data,
                              const UserInterface& UI){
    if (not(filename.empty()))
    {
        
        cout << "note: Writing output at " << filename ;
        if (filename == "ihixs_output") cout << "(the default output file)";
        cout << endl << endl;
        const char * output_fname = filename.c_str();
        fstream my_local_outfile(output_fname, fstream::out);
        if(my_local_outfile.is_open())
        {
            my_local_outfile.precision(5);
            my_local_outfile << "ihixs results " << endl;
            
            my_local_outfile << output << endl;
            my_local_outfile << endl << endl
            << output_text_separator() << endl;
            //
            my_local_outfile << "Results in Mathematica format" << endl;
            my_local_outfile<<endl<<data<<endl;
            //
            my_local_outfile << endl << endl
            << output_text_separator() << endl;
            //
            my_local_outfile << "User defined options (including command line ones)"<<endl;
            my_local_outfile << UI.GiveAllOptionsAndTheirValues();
            
            my_local_outfile << endl << endl << output_text_separator() <<endl;
            my_local_outfile << "Runcard used (options may have been overwritten by command line ones in runtime)"
                             << endl;
            my_local_outfile << "runcard_name=\""<<input_filename
            <<"\""<<endl;
            
            std::ifstream inputfile(input_filename, std::fstream::in);
            
            
        }
        else
        {
            cout<<"\nfailbit = "<<my_local_outfile.fail()<<endl;
            cout << "Error opening file "<<filename<<endl;
        }
        my_local_outfile.close();
    }
}
