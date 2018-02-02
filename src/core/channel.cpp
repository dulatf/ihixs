// Copyright 2015 ihixs team
#include<string>
#include<vector>

#include "src/core/channel.h"



ostream& operator<<(ostream& stream, const Channel& ch) {
    for (int i = 0; i < ch.NumberOfTerms(); i++)
        stream << *(ch._terms[i]);
    return stream;
}

Channel::Channel(const string& name) {
    _channel_name = name;
    checkName();
}

Channel::~Channel() {
    for (int i=0; i < _terms.size(); i++) {
        delete _terms[i];
    }
}


void Channel::AddTerm(const string& type,
                      const AsSeries& coeff,
                      LuminosityIntegral* LI) {
    const string term_name = _channel_name + " "+ type;
    _terms.push_back(new SigmaTerm(term_name, coeff, LI));
}



void Channel::checkName() {
    string implemented_channels[5] = {"gg", "qg", "qqbar", "qq", "q1q2"};
    bool found = false;
    for (int i = 0; i < 5; i++) {
        if (_channel_name == implemented_channels[i]) {
            found = true;
            break;
        }
    }
    if (not(found)) {
        cout << "\n[Error] you have requested an initial state channel "
            << "that is not implemented. Please try again!" << endl;
        exit(1);
    }
}

void Channel::Truncate(int int_qcd_perturbative_order) {
    for (int i = 0; i < _terms.size(); i++) {
        _terms[i]->Truncate(int_qcd_perturbative_order);
    }
}

void Channel::Evaluate(const InputParameters& input) {
    for (int i = 0; i < _terms.size(); i++) {
        _terms[i]->Evaluate(input);
    }
}

bool Channel::EvaluateSingleTerm(const InputParameters& input,
                                 const string& term_name) {
    bool found = false;
    for (int i = 0; i < _terms.size(); i++) {
        if (_terms[i]->Name() == term_name) {
            _terms[i]->Evaluate(input);
            found = true;
            break;
        }
    }
    return found;
}


string Channel::ListOfTerms() {
    stringstream stream;
    for (int i = 0; i < _terms.size(); i++) {
        stream << _terms[i]->Name() << endl;
    }
    return stream.str();
}

SigmaTerm* Channel::TermPtr(const string& name) {
    for (int i = 0; i < _terms.size(); i++) {
        if (_terms[i]->Name() == name) {
            return _terms[i];
        }
    }
    // not found
    cout << "Term " << name
        << " not found. Below is a list of terms in channel "
        << _channel_name << endl;
    cout << ListOfTerms() << endl;
    exit(0);
}



AsSeries Channel::QCDResult() {
    AsSeries res;
    for (int j = 0; j < _terms.size(); j++) {
        res = res + _terms[j]->QCDResult();
    }
    return res;
}

AsSeries Channel::EWResult() {
    AsSeries res;
    for (int j = 0; j < _terms.size(); j++) {
        res = res + _terms[j]->EWResult();
    }
    return res;
}

vector<double> Channel::ResultVector() {
    vector<vector<double> > result_vectors_from_sectors;
    for (int j = 0; j < _terms.size(); j++) {
        result_vectors_from_sectors.push_back(_terms[j]->ResultVector());
    }
    int number_of_components = result_vectors_from_sectors[0].size();
    vector<double> res(number_of_components, 0.0);
    for (int i = 0; i < number_of_components; i++) {
        for (int j = 0; j < _terms.size(); j++) {
            res[i] += result_vectors_from_sectors[j][i];
        }
    }
    return res;
}


vector<double> Channel::ErrorVector() {
    vector<vector<double> > result_vectors_from_sectors;
    for (int j = 0; j < _terms.size(); j++)  {
        result_vectors_from_sectors.push_back(_terms[j]->ErrorVector());
    }
    int number_of_components = result_vectors_from_sectors[0].size();
    vector<double> res(number_of_components, 0.0);
    for (int i = 0; i < number_of_components; i++) {
        for (int j = 0; j < _terms.size(); j++) {
            res[i] += pow(result_vectors_from_sectors[j][i], 2.);
        }
        res[i] = sqrt(res[i]);
    }
    return res;
}

