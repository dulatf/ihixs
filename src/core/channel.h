// Copyright 2015 ihixs team

#ifndef SRC_CORE_CHANNEL_H_
#define SRC_CORE_CHANNEL_H_

#include<string>
#include<vector>
#include "src/core/sigma_term.h"
#include "src/tools/as_series.h"



class Channel {
 public:
    explicit Channel(const string& name);
    ~Channel();
    void AddTerm(const string& type, const AsSeries& coeff,
                 LuminosityIntegral*);
    void Truncate(int m);
    void Evaluate(const InputParameters& input);
    bool EvaluateSingleTerm(const InputParameters& input, const string& name);
    AsSeries QCDResult();
    AsSeries EWResult();
    SigmaTerm* GiveTermPtr(int m) const {return _terms[m];}
    SigmaTerm* TermPtr(const string& name);
    int NumberOfTerms()const {return static_cast<int>(_terms.size());}
    string Name() const {return _channel_name;}
    friend ostream& operator<<(ostream&, const Channel&);
    string ListOfTerms();
    vector<double> ResultVector();
    vector<double> ErrorVector();

 protected:
    vector<SigmaTerm*> _terms;
    string _channel_name;

 private:
    void checkName();
};


#endif  // SRC_CORE_CHANNEL_H_

