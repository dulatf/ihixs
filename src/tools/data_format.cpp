#include "data_format.h"
#include<iomanip>// for setw
#include<iostream>// for cout

SingleDataEntry::SingleDataEntry() {
    _key = "";
    _res = "";
    _err = "";
    _val = ResultPair();
    _verbose_level = 2;
}

string SingleDataEntry::ConvertDoubleToString(const double & x){
    stringstream res_stream;
    res_stream << scientific;
    res_stream << setprecision(8);
    res_stream << x;
    return res_stream.str();
}
SingleDataEntry::SingleDataEntry(const string& key, const ResultPair& res, const int& verbose_level){
    _key=key;
    
    _res=ConvertDoubleToString(res.val());
    _err=ConvertDoubleToString(res.err());
    _val=res;
    _verbose_level = verbose_level;
}
SingleDataEntry::SingleDataEntry(const string& key, const string& res, const int& verbose_level){
    _key=key;
    _res=res;
    _err="";
    _verbose_level = verbose_level;
    _val = ResultPair();
    }
SingleDataEntry::SingleDataEntry(const string& key, const int& res, const int& verbose_level){
    _key=key;
    _res = std::to_string(res);
    _err="";
    _val=ResultPair(double(res),0.0);
    _verbose_level = verbose_level;
}
SingleDataEntry::SingleDataEntry(const string& key, const double& res, const int& verbose_level){
    _key=key;
    _res=ConvertDoubleToString(res);
    _err="";
    _val=ResultPair(res,0.0);
    _verbose_level = verbose_level;
}


string SingleDataEntry::MathematicaFormat(){
    stringstream datastring;
    datastring << scientific;
    datastring << setprecision(8);
    datastring<<"\""<<_key<<"\"->"<<_res;
    return datastring.str();
}

string SingleDataEntry::TextFormat(){
    stringstream datastring;
    
    if (!(_err==""))
        datastring<<setw(25)<<left<<_key<<" = "<<_res<<" ["<<_err<<"]";
    else
        datastring<<setw(25)<<left<<_key<<" = "<<_res;
    return datastring.str();
}




void DataFormat::Add(const string& key, const ResultPair& res){
    _entries.push_back(SingleDataEntry(key,res,0));
}

void DataFormat::Add(const string& key, const string& res){
    _entries.push_back(SingleDataEntry(key,res,0));
}

void DataFormat::Add(const string& key, const int& res){
    _entries.push_back(SingleDataEntry(key,res,0));
}

void DataFormat::Add(const string& key, const double& res){
    _entries.push_back(SingleDataEntry(key,res,0));
    
}

void DataFormat::Add(const string& key, const ResultPair& res, const string& verbose_level){
    _entries.push_back(SingleDataEntry(key,res,compute_integer_verbose_level(verbose_level)));
}

void DataFormat::Add(const string& key, const string& res, const string& verbose_level){
    _entries.push_back(SingleDataEntry(key,res,compute_integer_verbose_level(verbose_level)));
}

void DataFormat::Add(const string& key, const int& res, const string& verbose_level){
    _entries.push_back(SingleDataEntry(key,res,compute_integer_verbose_level(verbose_level)));
}

void DataFormat::Add(const string& key, const double& res, const string& verbose_level){
    _entries.push_back(SingleDataEntry(key,res,compute_integer_verbose_level(verbose_level)));
}

string DataFormat::str(){
    stringstream res;
    res<<"<|";
    for (int i=0;i<_entries.size()-1;i++){
        res<<_entries[i].MathematicaFormat()<<",";
    }
    res<<_entries[_entries.size()-1].MathematicaFormat();
    res<<"|>";
    return res.str();
}




string DataFormat::RestrictedOutput(){
    stringstream res;
    res << "Result" << endl;
    for (int i=0;i<_entries.size();i++){
        if ( _entries[i].VerboseLevel() < 1 )
            res<<_entries[i].TextFormat()<<endl;
    }
    return res.str();
}

string DataFormat::FullOutput(){
    stringstream res;
    res << "Result" << endl;
    for (int i=0;i<_entries.size();i++){
        if (_entries[i].VerboseLevel() < 2)
            res<<_entries[i].TextFormat()<<endl;
    }
    return res.str();
}

ResultPair DataFormat::GiveVal(const string& key) {
    for (int i = 0; i < _entries.size(); i++) {
        if (key == _entries[i].Key() ) {
            return _entries[i].Val();
        }
    }
    
    std::cout << "Error in DataFormat: you asked for an entry with a key: <" << key << "> that doesn't match any existing entry." << endl;
    exit(0);
}

bool DataFormat::Exists(const string& key) {
    for (int i = 0; i < _entries.size(); i++) {
        if (key == _entries[i].Key() ) {
            return true;
        }
    }
    return false;
    
    
}

int DataFormat::compute_integer_verbose_level(const string& verbose_level) {
    if (verbose_level == "silent") return 1;
    else if (verbose_level == "all") return 0;
    else if (verbose_level == "internal") return 2;
    else {
        cout << "Error in " << __func__ << " : verbosity level unknown <"
            << verbose_level << ">" << endl;
        exit(EXIT_FAILURE);
    }

}








