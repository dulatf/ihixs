

#include <iomanip>
#include <fstream>
#include <algorithm> // for remove_if in parser
#include <cctype> // for trim
#include <locale> // for trim
// C
#include <string.h> // for options of getopt
#include <getopt.h> // for options of getopt
#include <boost/algorithm/string/replace.hpp> // for string replacements

#include "user_interface.h"

int UIOption::get_type()
{
    if (type=="Required") return 1;
    else if (type=="Optional")  return 2;
    else if (type=="None") return 0;
    else
        {
        cout<<"\n wrong option"<<endl;
        exit(1);
        }
}



UserInterface::UserInterface()
{
    
    
    
    // operational options: a_s and perturbative order
    options.push_back( UIOption(
                                "Etot",
                               "COM energy of the collider in GeV",
                               "Required",
                               "13000.0",
                               "operational options"
                               ));
    //: masses and scales
    options.push_back(UIOption("m_higgs",
                               "higgs mass in GeV",
                               "Required",
                               "125.0",
                               "masses and scales options"));
    options.push_back( UIOption("mur",
                                "mur",
                                "Required",
                                "62.5",
                                "masses and scales options"));
    options.push_back( UIOption("muf",
                                "muf",
                                "Required",
                                "62.5",
                                "masses and scales options"));
    
    
    
    
    //: operational options: PDFs
    
    options.push_back( UIOption("pdf_member",
                                "pdf member id (the range depends on the pdf set)",
                                "Required",
                                "0",
                                "operational options"));
    
    options.push_back( UIOption("pdf_set",
                                "choose a specific pdf set name (LHAPDF6 list at lhapdf.hepforge.org/pdfsets.html). This set will be used irrespectively of order.",
                                "Required",
                                "PDF4LHC15_nnlo_100",
                                "operational options"));
    
    options.push_back( UIOption("pdf_set_for_nlo",
                                "pdf set used when computing PDF-TH error.",
                                "Required",
                                "PDF4LHC15_nlo_100",
                                "operational options"));
    
    
    
    //: operational options: what to compute
    
    options.push_back( UIOption("with_eft",
                                "compute the cross section in the EFT approximation",
                                "Required",
                                "true",
                                "operational options"));
    
    
    
    options.push_back( UIOption("with_exact_qcd_corrections",
                                "true to include the exact quark mass effects at NLO, false to omit them",
                                "Required",
                                "false",
                                "operational options"));
    
    options.push_back( UIOption("with_ew_corrections",
                                "true to include the exact quark mass effects at NLO, false to omit them",
                                "Required",
                                "false",
                                "operational options"));
    
    options.push_back( UIOption(
                                "with_mt_expansion",
                                "include NNLO 1/mt terms",
                                "Required",
                                "false",
                                "operational options"));
    
    
//    options.push_back( UIOption(
//                                "with_truncation_error",
//                                "compute truncation error",
//                                "Required",
//                                "false",
//                                "operational options"));
    options.push_back( UIOption(
                                "with_delta_pdf_th",
                                "compute PDF-TH uncertainty",
                                "Required",
                                "false",
                                "operational options"));
    
    options.push_back( UIOption("with_scale_variation",
                                "estimate scale variation (mur and muf should be at mh/2)",
                                "Required",
                                "false",
                                "operational options"));
    options.push_back( UIOption("with_indiv_mass_effects",
                                "compute separately light quark contributions",
                                "Required",
                                "false",
                                "operational options"));
    options.push_back( UIOption("with_pdf_error",
                                "whether or not to compute error due to pdfs",
                                "Required",
                                "false",
                                "operational options"));
    options.push_back( UIOption(
                                "with_a_s_error",
                                "compute a_s uncertainty",
                                "Required",
                                "false",
                                "operational options"));
    options.push_back( UIOption(
                                "with_resummation",
                                "include threshold resummation",
                                "Required",
                                "false",
                                "operational options"));
    
    
    
    // operational options: resummation
    options.push_back( UIOption("resummation_log_order",
                                "0:LL, 1:NLL, 2:NNLL, 3:N3LL",
                                "Required",
                                "3",
                                "operational options"));// default N3LL
    options.push_back( UIOption("resummation_matching_order",
                                "0:L0, 1:NL0, 2:NNL0, 3:N3L0",
                                "Required",
                                "3",
                                "operational options"));// default N3LO
    options.push_back( UIOption("resummation_type",
                                "variant of threshold resummation, i.e. log:classical, psi, AP2log, AP2psi ",
                                "Required",
                                "log",
                                "operational options"));
    
    options.push_back( UIOption(
                                "with_scet",
                                "include scet resummation",
                                "Required",
                                "false",
                                "operational options"));
    
    // IO options
    options.push_back( UIOption("verbose",
                                "level of verbosity: minimal or medium. Medium shows channel breakdown EFT cross section.",
                                "Required",
                                "minimal",
                                "input-output options",
                                'v'));
    options.push_back( UIOption("input_filename",
                                "filename to use as runcard",
                                "Required",
                                "default.card",
                                "input-output options",
                                'i'));
    options.push_back( UIOption("output_filename",
                                "filename to write output",
                                "Required",
                                "ihixs_output",
                                "input-output options",
                                'o'));
    options.push_back( UIOption("help",
                                "print all options and help messages per option.",
                                "Optional",
                                "false",
                                "input-output options",
                                'h'));
    options.push_back( UIOption("make_runcard",
                                "create default runcard file as default_card.",
                                "Optional",
                                "false",
                                "input-output options",
                                'd'));
    options.push_back( UIOption("make_pheno_card",
                                "create pheno runcard file as pheno_card.",
                                "Optional",
                                "false",
                                "input-output options",
                                'p'));
    options.push_back( UIOption("write_documentation",
                                "print the help message in a TeX form.",
                                "Optional",
                                "false",
                                "input-output options",
                                'w'));
    options.push_back( UIOption("with_eft_channel_info",
                                "print eft cross section per channel per order",
                                "Optional",
                                "false",
                                "input-output options"));
    options.push_back( UIOption("with_resummation_info",
                                "info from resummation: true, false",
                                "Required",
                                "false",
                                "input-output options"));
    
    //: masses and scales
    //: m_t(m_t) = 163.7 is the pdg equivalent
    //: m_t(m_t) = 162.7 default based on HXSWG internal note (to be compatible with ATLAS & CMS values)
    options.push_back( UIOption("mt_msbar",
                                "MSbar top mass",
                                "Required",
                                "162.7",
                                "masses and scales options"));
    options.push_back( UIOption("mt_msbar_ref_scale",
                                "reference scale for the top mass in MSbar",
                                "Required",
                                "162.7",
                                "masses and scales options"));
    options.push_back( UIOption("mt_on_shell",
                                "On Shell top mass",
                                "Required",
                                "172.5",
                                "masses and scales options"));
    //: m_b defaults following http://arxiv.org/pdf/0907.2110.pdf
    options.push_back( UIOption("mb_msbar",
                                "MSbar bottom mass",
                                "Required",
                                "4.18",
                                "masses and scales options"));
    options.push_back( UIOption("mb_msbar_ref_scale",
                                "reference scale for the bottom mass in MSbar",
                                "Required",
                                "4.18",
                                "masses and scales options"));
    //: m_b on-shell HXSWG internal note (using the four loop conversion formula from the MSBAR mass)
    options.push_back( UIOption("mb_on_shell",
                                "On Shell bottom mass",
                                "Required",
                                "4.92",
                                "masses and scales options"));

    //: m_c defaults following http://arxiv.org/pdf/0907.2110.pdf
    options.push_back( UIOption("mc_msbar",
                                "MSbar charm mass",
                                "Required",
                                "0.986",
                                "masses and scales options"));
    options.push_back( UIOption("mc_msbar_ref_scale",
                                "reference scale for the charm mass in MSbar",
                                "Required",
                                "3.0",
                                "masses and scales options"));
    //:m_c on-shell from pdg http://pdg.lbl.gov/2014/listings/rpp2014-list-c-quark.pdf
    //: also HXSWG internal note (from the mb-mc formula hep-ph/0408002)
    options.push_back( UIOption("mc_on_shell",
                                "On Shell charm mass",
                                "Required",
                                "1.67",
                                "masses and scales options"));

    options.push_back( UIOption("top_scheme",
                                "msbar or on-shell",
                                "Required",
                                "msbar",
                                "masses and scales options"));
    options.push_back( UIOption("bottom_scheme",
                                "msbar or on-shell",
                                "Required",
                                "msbar",
                                "masses and scales options"));
    options.push_back( UIOption("charm_scheme",
                                "msbar or on-shell",
                                "Required",
                                "msbar",
                                "masses and scales options"));

    //: Yukawas relative to the SM one
    options.push_back( UIOption("y_top",
                                "factor multiplying the Yt. Set to zero to remove the top quark",
                                "Required",
                                "1.0",
                                "masses and scales options"));
    options.push_back( UIOption("y_bot",
                                "factor multiplying the Yb. Set to zero to remove the bottom quark",
                                "Required",
                                "1.0",
                                "masses and scales options"));
    options.push_back( UIOption("y_charm",
                                "factor multiplying the Yc. Set to zero to remove the charm quark",
                                "Required",
                                "1.0",
                                "masses and scales options"));
    
    //: widths
    options.push_back( UIOption("gamma_top",
                                "width of top quark",
                                "Required",
                                "0.0",
                                "masses and scales options"));

    options.push_back( UIOption("gamma_bot",
                                "width of bottom quark",
                                "Required",
                                "0.0",
                                "masses and scales options"));

    options.push_back( UIOption("gamma_charm",
                                "width of charm quark",
                                "Required",
                                "0.0",
                                "masses and scales options"));
    
    // operational options: a_s and perturbative order
    options.push_back( UIOption("qcd_perturbative_order",
                                "LO, NLO, NNLO, N3LO : ihixs will compute up to this order in a_s",
                                "Required",
                                "N3LO",
                                "operational options"));
    options.push_back( UIOption("with_fixed_as_at_mz",
                                "set the value of a_s(mZ) by hand. Beware: this might not be compatible with your pdf choice.",
                                "Required",
                                "0.0",
                                "operational options"));
    options.push_back( UIOption("qcd_order_evol",
                                "used for a_s and quark mass evolution 0:L0, 1:NL0, 2:NNL0, 3:N3L0",
                                "Required",
                                "3",
                                "operational options"));// default N3LO
    options.push_back( UIOption(
                                "with_lower_ord_scale_var",
                                "also compute scale variation for lower than the current order",
                                "Required",
                                "false",
                                "operational options"));
    //: numerical precision options
    
    options.push_back( UIOption("epsrel",
                                "cuba argument: target relative error",
                                "Required",
                                "0.0001",
                                "numerical precision options"));
    options.push_back( UIOption("epsabs",
                                "cuba argument: target absolute error",
                                "Required",
                                "0.0",
                                "numerical precision options"));
    options.push_back( UIOption("mineval",
                                "cuba argument: minimum points to be evaluated",
                                "Required",
                                "50000",
                                "numerical precision options"));
    options.push_back( UIOption("maxeval",
                                "cuba argument: maximum points to be evaluated",
                                "Required",
                                "50000000",
                                "numerical precision options"));
    options.push_back( UIOption("nstart",
                                "cuba argument: number of points for first iteration",
                                "Required",
                                "10000",
                                "numerical precision options"));
    options.push_back( UIOption("nincrease",
                                "cuba argument: number of points for step increase",
                                "Required",
                                "1000",
                                "numerical precision options"));
    options.push_back( UIOption("cuba_verbose",
                                "cuba argument: verbosity level: 0=silent, 2=iterations printed out",
                                "Required",
                                "0",
                                "numerical precision options"));

    
    
    
    
}

void UserInterface::SetOption(const string& key, const string& val) {
    for (int i = 0 ; i < options.size(); i++) {
        if (options[i].name == key) {
            options[i].value = val;
            break;
            }
    }
}

void UserInterface::SetOption(const string& key, const double& val) {
    stringstream s;
    s<<val;
    SetOption(key,s.str());
}

void UserInterface::ParseInput(int argc, char * const *argv)
{
     // parse command line arguments
     vector<vector<string> > parsed_options = ParseCmd(argc, argv, true);
     // look for runcard definition in command line (-i <cardname>)
     // or for --help
     for (unsigned i=0;i<parsed_options.size();i++)
          {
          if (parsed_options[i][0]=="input_filename")
               {
               SetOption("input_filename",parsed_options[i][1]);
               }
          if (parsed_options[i][0]=="help")
              {
                  print_help_message();
                  exit(0);
              }
          if (parsed_options[i][0]=="make_runcard")
              {
                  MakeRuncard("default.card");
                  exit(0);
              }
              if (parsed_options[i][0]=="make_pheno_card")
              {
                  MakePhenoCard();
                  exit(0);
              }
        if (parsed_options[i][0]=="write_documentation") {
                  WriteDocumentation();
                  exit(0);
              }

          }
     // Read the (potentially user defined) runcard
     ParseFile(giveString("input_filename"), true);
    // Modify parameters that were declared in command line,
    // overwritting those of runcard
     for (unsigned i=0;i<parsed_options.size();i++)
          {
          for (unsigned j=0;j<options.size();j++)
               {
               if (parsed_options[i][0]==options[j].name)
                    {
                        options[j].value=parsed_options[i][1];
                    break;
                    }
               }

          }
    //cout << "***\n" << GiveAllOptionsAndTheirValues() << endl;exit(0);
}


string UserInterface::giveString(const string& option_name) const {
    for (int i=0; i < options.size(); i++) {
        if (options[i].name == option_name) return options[i].value;
    }
    cout << "Error in UserInterface : option not found : " << option_name <<endl;
    exit(1);
}

double UserInterface::giveDouble(const string& option_name) const {
    for (int i=0; i < options.size(); i++) {
        if (options[i].name == option_name) return atof(options[i].value.c_str());
    }
    cout << "Error in UserInterface : option not found : " << option_name <<endl;
    exit(1);
}

int UserInterface::giveInt(const string& option_name) const {
    for (int i=0; i < options.size(); i++) {
        if (options[i].name == option_name) return atoi(options[i].value.c_str());
    }
    cout << "Error in UserInterface : option not found : " << option_name <<endl;
    exit(1);
}

bool UserInterface::giveBool(const string& option_name) const {
    for (int i=0; i < options.size(); i++) {
        if (options[i].name == option_name) {
            if (options[i].value == "true") return true;
            else if (options[i].value == "True") return true;
            else if (options[i].value == "false") return false;
            else if (options[i].value == "False") return false;
            else {
                cout << "Error in UserInterface: variable "
                    << option_name
                    << "is a Boolean and was set to the unrecognized "
                    << options[i].value << endl;
                exit(1);
            }
        }
    }
    cout << "Error in UserInterface : option not found : " << option_name <<endl;
    exit(1);
}


int UserInterface::ParseFile(const string& in, bool verbose)
{
    // oooopen
    std::fstream file(in.c_str(), std::fstream::in);


    // and do the job
    if(!file.good()) {
        stringstream ss;
        if (in == "default.card") {
        std::cout << "Tried to open the default runcard: " << in << " but failed!" << endl;
        cout << "The program will use default parameters overwritten by command line options." << std::endl;
        cout << "The parameter values used can be found in the output file (see below for the filename)" << endl;
        return 1;
        }
        else {
            std::cout << "Tried to open the runcard you specified: " << in << " but failed!" << endl;
            std::cout << "Please make sure that " << in << " is in the current directory." << endl;
            exit(1);
            
        }
    }
    else{
        char buff[256];
        // get first line
        file.getline(buff, 256);
        // loop over lines
        for(unsigned c=1; !file.eof(); ++c){
            // convert to string
            std::string s;
            // ...discard spaces and comments
            for(char *p=buff; *p!='\0'; ++p){
                if(*p == '#')break;
                else{
                    if(*p != ' ') s += *p;
                }
            }
            // removing empty spaces again
            s.erase( remove_if( s.begin(),s.end(),::isspace), s.end() );
            
            // the line below is reported to be more generic than ::isspace
            //std::bind(std::isspace<char>, _1, std::locale::classic() )
            // non-empty lines
            if(!s.empty()){
                //cout<<"<"<<s<<">"<<endl;

                size_t pos = s.find('=');
                // Invalid line ( a line without an `=` sign)
                if(pos == s.npos){
                    if(verbose) {
                        std::cout << "Invalid line in input file " << in << "line :" << c << "<" << s <<">"<< std::endl;
                        exit(1);
                    }
                }
                else{
                    // valid line !
                    bool did=false;
                    for(unsigned i=0; i<options.size(); ++i)
                        if(s.substr(0, pos) == options[i].name){
                            string value =s.substr(pos+1, s.npos);
                            //cout<<"value=<"<<value<<">"<<endl;
                            if (!value.empty()){
                                //options[i].set(value);
                                //new
                                options[i].value=value;
                            }
                            else{
                                cout<<"Option "<<options[i].name<<" recognized but it is set to no value."
                                    <<endl;
                                exit(1);
                            }

                            did=true;
                            break;
                        }

                    if(!did && verbose){
                        stringstream msg;
                        msg << "Unrecognised option in input file " << in << ": " << s.substr(0,pos);
                        throw(msg.str());
                    }
                }
            }

            // iterate
            file.getline(buff, 256);
        }
    }
     return 0;
}

option * UserInterface::create_getopt_option_array()
{
     // create the long_option array
     // it has #options elements
     int N=int(options.size());
     option *long_options = new option[N+1];
     // feed them
     for(unsigned i=0; i<N; ++i)
          {
          long_options[i].name = new char[strlen(options[i].name.c_str())];
          strcpy(const_cast<char*>(long_options[i].name), options[i].name.c_str());
          long_options[i].has_arg = options[i].get_type();
          long_options[i].flag = NULL;
          long_options[i].val = options[i].short_name;
          }
     // last element (viva old school C)
     long_options[N+1].name = NULL;
     long_options[N+1].has_arg = 0;
     long_options[N+1].flag = NULL;
     long_options[N+1].val = 0;


     return long_options;
}

string UserInterface::create_getopt_optdesc()
{
     // create short options descriptor
     std::string optdesc;
     for(unsigned i=0; i<options.size(); ++i)
          {
          // help and version need special treatment later
          if(options[i].short_name!=0) //: 0 is the default value for short_name when there is no short name
               {
               // feed the optiondesc string
               optdesc += options[i].short_name;
               if(options[i].get_type() == required_argument) optdesc += ":";
               if(options[i].get_type() == optional_argument) optdesc += "::";
               }
          }
     return optdesc;
}

void UserInterface::print_help_message()
{
    std::cout
    << "Input parameters and settings (\"options\") in ihixs are defined either "
    << std::endl
    << "in the runcard or from command line. Command line options have to be"
    << std::endl
    << "given the gnu way, i.e. '--<option> <value>' or '-<shorthand> <value>'. "
    << std::endl
    << "For example, to set the collider energy to 14TeV use:"
    << std::endl << std::endl
    << "\tihixs  --Etot 14000 "
    << std::endl << std::endl
    << "Shorthands exist for selected options, and are denoted below by [<shorthand>]. " << std::endl
    << "For example to run ihixs with input card my.card and output file my.out"
    << std::endl << std::endl
    << "\tihixs -i my.card -o my.out" << std::endl << std::endl
    << "Please note that --option = value doesn't work!!"
    << std::endl << std::endl
    << "Here is a list of available options: "
    << std::endl << std::endl;
     // Print all other options (version is included here)
     for(unsigned i=0; i<options.size(); i++)
          {
          stringstream option_name;
          option_name << options[i].name;
          if ( options[i].short_name != 0 )
              option_name << " [" << options[i].short_name << "]";
          std::cout << ' ' << std::setw(25) << std::left
          << option_name.str() << " "<< options[i].desc << std::endl;
          }
     exit(1);
}


string UserInterface::GiveAllOptionsAndTheirValues() const {
    stringstream res;
    for(unsigned i=0; i<options.size(); i++)
    {
        
        res << std::setw(27) << std::left
            << options[i].name
            << " = "
            << std::setw(22) << std::left
            << options[i].value
            << std::setw(4) << std::left
            << "#"
            //<< std::setw(50) << std::left
            << " " << options[i].desc
            << std::endl;
    }
    return res.str();
}

void UserInterface::MakeRuncard(const char * output_fname)
{
    stringstream runcard_text;
    runcard_text << "# " << output_fname  << endl
    << "# Command-line options can override the runcard options below."
    << std::endl
    << "# The final values used by the program can be seen at the output file."
    << std::endl
    << "# Anything after the character '#' in a line is considered a comment and is ignored by the parser." << endl
    << "# In case of having unwillingly modify this card, you can regenerate it by ./ihixs -d"
    << std::endl << std::endl;
    // Print all other options (version is included here)
    runcard_text << GiveAllOptionsAndTheirValues();
    //const char * output_fname = "default_card";
    fstream my_local_outfile(output_fname, fstream::out);
    if(my_local_outfile.is_open()){
        my_local_outfile << runcard_text.str();
        cout << "Creating  : " << output_fname
             << " and exiting. " << endl;
    }
    else
    {
        cout << "\nfailbit = " << my_local_outfile.fail() << endl;
        cout << "Error opening file " << output_fname << endl;
    }
    my_local_outfile.close();
    exit(1);
}

void UserInterface::MakePhenoCard() {
    SetOption("with_exact_qcd_corrections","true");
    SetOption("with_ew_corrections","true");
    SetOption("with_mt_expansion","true");
    SetOption("with_delta_pdf_th","true");
    SetOption("with_scale_variation","true");
    SetOption("with_pdf_error","true");
    SetOption("with_a_s_error","true");
    MakeRuncard("pheno.card");
}


void UserInterface::WriteDocumentation()
{
    stringstream tex_output;
    tex_output << "\\begin{table}" << endl
    << "\\bgfb" << endl
    << "\\multicolumn{2}{c}{\\textbf{Table~\\ref{tab:operationaloptions}: Operational options}}" << endl
    << "\\\\ "
    << endl;
    tex_output << WriteTexDocForClassification("operational options");
    tex_output << "\\egfb" << endl
    << "\\textcolor{white}{\\caption{\\label{tab:operationaloptions}}}" << endl
    << "\\end{table}" << endl;
    
    tex_output << "\\begin{table}" << endl
    << "\\bgfb" << endl
    << "\\multicolumn{2}{c}{\\textbf{Table~\\ref{tab:inputoutputoptions}: Input-Output options}}" << endl
    << "\\\\ "
    << endl;
    tex_output << WriteTexDocForClassification("input-output options");
    tex_output << "\\egfb" << endl
    << "\\textcolor{white}{\\caption{\\label{tab:inputoutputoptions}}}" << endl
    << "\\end{table}" << endl;
    
    tex_output << "\\begin{table}" << endl
    << "\\bgfb" << endl
    << "\\multicolumn{2}{c}{\\textbf{Table~\\ref{tab:massesandscalesoptions}: Masses and scales options}}" << endl
    << "\\\\ "
    << endl;
    tex_output << WriteTexDocForClassification("masses and scales options");
    tex_output << "\\egfb" << endl
    << "\\textcolor{white}{\\caption{\\label{tab:massesandscalesoptions}}}" << endl
    << "\\end{table}" << endl;
    
    tex_output << "\\begin{table}" << endl
    << "\\bgfb" << endl
    << "\\multicolumn{2}{c}{\\textbf{Table~\\ref{tab:numericalprecisionoptions}: Numerical precision options}}" << endl
    << "\\\\ "
    << endl;
    tex_output << WriteTexDocForClassification("numerical precision options");
    tex_output << "\\egfb" << endl
    << "\\textcolor{white}{\\caption{\\label{tab:numericalprecisionoptions}}}" << endl
    << "\\end{table}" << endl;
    
    string TeX_string=tex_output.str();
    boost::replace_all(TeX_string, "_", "\\_");
    cout << TeX_string << endl;
    exit(1);
}

string UserInterface::WriteTexDocForClassification(const string& classification_name) {
    stringstream tex_output;
    for(unsigned i=0; i<options.size(); i++) {
        if (options[i].classification==classification_name) {
            tex_output << "\\tt{"
            << options[i].name << "} : \\textbf{"
            << options[i].value
            << "} & "
            << options[i].desc
            << std::endl
            << "\\\\"
            << std::endl;
        }
    }
    return tex_output.str();
}


void UserInterface::PrintAllOptions() const
{
    cout<<"\n-----------------------------------------------------------------";
    cout<<"UI options:"<<endl;
    for (int i=0;i<options.size();i++)
        {
        cout<<"\n"<<options[i].name<<" : "<<options[i].print();
        }
    cout<<"\n-----------------------------------------------------------------";
    cout<<endl;
}

string UserInterface::PrintOptionsAndValues() const
{
    stringstream st;
    for (int i = 0; i < options.size(); i++ ) {
        st << options[i].name << " : " << options[i].print() << endl;
    }
    return st.str();
}

// trim from start (in place)
static inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
        return !std::isspace(ch);
    }));
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}


vector<vector<string> > UserInterface::ParseCmd(int argc,  char * const *argv, bool locverbose)
{
    //cout<<"command line args: ";
    //for (int i=0;i<argc;i++) cout<<argv[i]<<" ";
    //cout<<endl;
     // getopt quirckiness:
     char * const * argv_for_getopt = argv;
     // create the long_option array
     option *long_options = create_getopt_option_array();
     // create short options descriptor
     string optdesc=create_getopt_optdesc();
     // we will put all found options in parsed_options
     vector<vector<string> > parsed_options;

     // Get all options with getopt_long
     int option_index=0;
     int c;
     // reinit global (?) vars... (i do not get it but seems to work)

     optarg=0;
     optind=0;
     //: getopt_long returns a character equal to
     //: 0 for any other LONG option
     //: -1 for the end of parsing
     //: 'x' for any short option corresponding to short-name='x'
     //: in the case 0, the variable option_index holds
     //: the position of the recognized option in the long_option array
     while((c = getopt_long(argc, argv_for_getopt, optdesc.c_str(), long_options, &option_index)) != -1)
          {
          // do flag
          //bool did=false;
          switch (c)
               {
                    // --long_option
                    case 0:
                    // loop over other possibilities
                    for(unsigned i=0; i<options.size(); ++i)
                         {
                         // check if indeed the name in long_options is the same as the name in options
                         if( options[i].name == long_options[option_index].name )
                              {
                              vector<string> loc_parsed_opt;
                              loc_parsed_opt.push_back(options[i].name);
                              //did = true;


                              // for normal options
                              if(optarg!=0)
                                   {
                                   string locarg = optarg;
                                   loc_parsed_opt.push_back(optarg);
                                   }
                              // for flags
                              else loc_parsed_opt.push_back("true");//options[i].set("1");
                              parsed_options.push_back(loc_parsed_opt);
                              break;
                              }
                         }

                    break;

                    case '?':// getopt_long already printed an error message.
                       {stringstream msg;
                       msg << "Unrecognised command line option";
                       throw(msg.str());

                    //exit(1);
                    break;}
                    // short options
                    default:
                    // loop over expected short options
                    for(unsigned i=0; i<options.size(); ++i)
                         {
                         if( c == options[i].short_name )
                              {
                              vector<string> loc_parsed_opt;
                              loc_parsed_opt.push_back(options[i].name);
                              //did = true;
                              // for normal options
                              if(optarg)
                                   {
                                   loc_parsed_opt.push_back(optarg);
                                   }
                              // for flags
                              else loc_parsed_opt.push_back("true");
                              parsed_options.push_back(loc_parsed_opt);
                              break;
                              }
                         }
                    }
          }

    for (int i=0; i < parsed_options.size(); i++) {
        if (parsed_options[i][1].find("=") != std::string::npos) {
            cout << "Error in parsing : you probably used the format --" << parsed_options[i][0]
            << " = X which doesn't work." << endl
            << "Please try --" << parsed_options[i][0] << " X (i.e. with an empty space instead of '=')" << endl
            << "The format " << parsed_options[i][0] << "=X also works (i.e. without empty spaces around '=')."
            << endl;
            exit(0);
        }
        trim(parsed_options[i][1]);
        if (parsed_options[i][1].empty()) {
            cout << "Error in parsing : you probably used the format --" << parsed_options[i][0]
            << "= X which doesn't work." << endl
            << " Please try --" << parsed_options[i][0] << " X (i.e. with an empty space instead of '=')" << endl;
            cout << "You might also have forgotten to define " << parsed_options[i][0]
                << " which is also illegal. "<< endl
                << "The format " << parsed_options[i][0] << "=X also works (i.e. without empty spaces around '=')."
            << endl;

            exit(0);
        }
    }
     return parsed_options;
}
