#ifndef GGF_MANAGER_H
#define GGF_MANAGER_H

#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <vector>
#include <iostream>
#include <iomanip>
#include <map>
#include <exception>
#include <math.h>
#include "inclusive_process.h"
#include "inclusive_higgs_eft.h"
#include "inclusive_higgs_eft_fast.h"
#include "inclusive_mt_expansion.h"
#include "inclusive_higgs_exact.h"
#include "threshold_resummation.h"
#include "data_format.h"// for DataFormat (mathematica data structure output)
#include "user_interface.h"
using namespace std;


#include "LHAPDF/LHAPDF.h"
#include "threshold_resummation.h"

#include "psi.h"


#include "scetter.h"

class GgfManager{
public:
    GgfManager(const UserInterface& UI);
    
private:
    DataFormat mydata;
    stringstream output;
    stringstream data;
    stringstream messages;
    
    InputParametersForThresRes* _tr_input;
    ResultPair _Higgs_XS;
    int _qcd_perturbative_order_int;
    bool _we_do_mt_expansion;
    bool _perform_pdf_error;
private:
    void perform_eft(const UserInterface& UI);

    void perform_mt_expansion(const UserInterface& UI);
    void perform_qcd_exact(const UserInterface& UI);
    void perform_qcd_exact_study(const UserInterface& UI);
    void perform_qcd_exact_plain(const UserInterface& UI);
    
    void perform_scet(const UserInterface& UI);
    
    void perform_catani_resummation(const UserInterface& UI);
    void perform_pdf_th_error(const UserInterface& UI);
    void perform_as_error(const UserInterface& UI);

    void perform_exact_nlo_with_selected_quarks(const UserInterface& UI, double yt, double yb, double yc);
    
    void compute_total_xs(const UserInterface& UI);
    void compute_uncertainty(const UserInterface& UI);
    ResultPair compute_delta_tbc(const UserInterface& UI);
    ResultPair compute_delta_ew(const UserInterface& UI);
    ResultPair compute_delta_trunc(const UserInterface& UI);
    ResultPair compute_delta_pdf_th(const UserInterface& UI);


    void single(UserInterface UI);
    void perform_pdf_error(UserInterface UI);
    void truncation_error(const UserInterface& UI);
    void full_logs_vs_truncated_logs_error(const UserInterface& UI);
    
    void write_output(const string& content,
                      const string& output_fname,
                      const string& input_fname,
                      const string& data,
                      const UserInterface& UI);
    
    int determine_int_qcd_perturbative_order(const string& );
    vector<double> compute_scale_variation_boundary(const string order);
    void perform_scale_variation(const UserInterface& UI);
    void compute_scale_variation_for_unrescaled_eft(InclusiveHiggsEFT* low, InclusiveHiggsEFT* high,int);
    string make_string_for_scale_variation_boundary(const double& val);
    string output_text_separator();
    
};

#endif

