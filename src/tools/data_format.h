#ifndef DATA_FORMAT_H
#define DATA_FORMAT_H


#include <sstream> // for stringstream

#include <vector>
#include "result_pair.h"
using namespace std;
// data formatted for Mathematica plots


class SingleDataEntry{
public:
    SingleDataEntry();
    SingleDataEntry(const string& key, const ResultPair& res, const int& verbose_level);
    SingleDataEntry(const string& key, const string& res, const int& verbose_level);
    SingleDataEntry(const string& key, const int& res, const int& verbose_level);
    SingleDataEntry(const string& key, const double& res, const int& verbose_level);
    // printing functions
    string MathematicaFormat();
    string TextFormat();
    // 
    ResultPair Val(){return _val;}
    string Key(){return _key;}
    int VerboseLevel() {return _verbose_level;}
private:
    string _key;
    string _res;
    string _err;
    ResultPair _val;
    int _verbose_level;
};

class DataFormat{
public:
    void Add(const string& key, const ResultPair& val);
    void Add(const string& key, const string& val);
    void Add(const string& key, const int& val);
    void Add(const string& key, const double& val);
    
    void Add(const string& key, const ResultPair& val, const string& verbose_level);
    void Add(const string& key, const string& val, const string& verbose_level);
    void Add(const string& key, const int& val, const string& verbose_level);
    void Add(const string& key, const double& val, const string& verbose_level);
    
    ResultPair GiveVal(const string& key);
    bool Exists(const string& key);
    string str();
    string RestrictedOutput();
    string FullOutput();
    //string Output(const string& verbose_level);
private:
    vector<SingleDataEntry> _entries;
    int compute_integer_verbose_level(const string& verbose_level);
    
    
};
#endif
