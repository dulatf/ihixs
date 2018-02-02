#ifndef USER_INTERFACE_H
#define USER_INTERFACE_H

#include<string>
#include<vector>
#include<iostream>
#include <stdlib.h> //: for exit()
#include <sstream> //: for stringstream

using namespace std;

//: declaring the getopt specialized class option (note small initial 'o')
struct option;

/** UIOption descriptor.
 *
 * Encloses all the information relevant for the UIOption parsing. */
class UIOption {
public:
    UIOption(const std::string& name_,
            const std::string& desc_,
            const string& type_,
            const string& default_,
            const string& classification_,
            char short_name_)
            :   type(type_),
                name(name_),
                short_name(short_name_),
                desc(desc_),
                value(default_),
                classification(classification_){};

    UIOption(const std::string& name_,
             const std::string& desc_,
             const string& type_,
             const string& default_,
             const string& classification_)
            :   type(type_),
            name(name_),
            short_name(0),
            desc(desc_),
            value(default_),
            classification(classification_){};
    ~UIOption(){};
    int get_type();
    /** Set  */
    void set(const string& val) {value = val;}
    /** Print... */
    string print() const {return value;}
    
public:
    string type;
    std::string name;
    char short_name;
    std::string desc;
    std::string value;
    std::string classification;
};




class UserInterface
{
public://methods
    UserInterface();
    ~UserInterface(){};
    void ParseInput(int argc, char * const *argv);
    void print_help_message();
    
    void PrintAllOptions() const;
    string PrintOptionsAndValues() const;
    string GiveAllOptionsAndTheirValues() const;
    
    string giveString(const string& option_name) const;
    double giveDouble(const string& option_name) const;
    int giveInt(const string& option_name) const;
    bool giveBool(const string& option_name) const;

    void SetOption(const string& key, const string& val);
    void SetOption(const string& key, const double& val);

public:
private://methods
    int ParseFile(const string &, bool);
    vector<vector<string> > ParseCmd(int argc,  char* const *argv,
                                        bool verbose);
     //: getopt interface
     option * create_getopt_option_array();
     string create_getopt_optdesc();
    
    void MakePhenoCard();
    void MakeRuncard(const char * output_fname);
    void WriteDocumentation();
    string WriteTexDocForClassification(const string& classification_name);

private://data
    vector<UIOption> options;
};

#endif
