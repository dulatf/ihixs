// Copyright 2015 ihixs team

#include<string>
#include "src/core/inclusive_process.h"
#include "src/tools/constants.h"



ostream& operator<<(ostream& stream, const InclusiveProcess& ip) {
    for (int i = 0; i < ip._channels.size(); i++) {
        stream << *(ip._channels[i]);
    }
    return stream;
}


InclusiveProcess::InclusiveProcess(const UserInterface& UI) {
    Initialize(UI);
}

void InclusiveProcess::Initialize(const UserInterface& UI) {
    _input.Configure(UI);
}

InclusiveProcess::~InclusiveProcess() {
    for (int i=0; i < _channels.size(); i++)
        delete _channels[i];
}

void InclusiveProcess::Evaluate() {
    this->SetUpContributions();
    for (int i = 0; i < _channels.size(); i++) {
        cout << "computing " << _channels[i]->Name()
            << " channel " << endl;
        _channels[i]->Evaluate(_input);
    }
    this->EvaluatePostVegasDependingOnProcess();
}

void InclusiveProcess::EvaluateSingleTerm(const string& term_name) {
    bool found = false;
    for (int i = 0; i < _channels.size(); i++) {
        if ( _channels[i]->EvaluateSingleTerm(_input, term_name) ) {
            found = true;
            break;
        }
    }
    this->EvaluatePostVegasDependingOnProcess();
    if ( !(found) ) {
        cout << "In EvaluateSingleTerm: no \"" << term_name
                     << "\" found in any  channel " << endl;
        exit(0);
    }
}



string InclusiveProcess::ChannelBreakdown() {
    stringstream st;
    st << "------------------------" << endl;
    st << "XS per channel" << endl;
    for (int i = 0; i < _channels.size(); i++) {
        AsSeries res = _channels[i]->QCDResult();
        st << setw(6) << left << _channels[i]->Name() << ":"
            << setw(20) << res.AddUp() << " = "
            << setw(20) << res.term_of_order(2)
            << setw(20) << res.term_of_order(3)
            << setw(20) << res.term_of_order(4)
            << setw(20) << res.term_of_order(5)
            << endl;
    }
    return st.str();
}


Channel* InclusiveProcess::GiveChannel(const string& name) {
    for (int i = 0; i < _channels.size(); i++) {
        if (_channels[i]->Name() == name) {
            return _channels[i];
        }
    }
    // not found
    cout << "Channel " << name
        << " not found. Below is a list of channels defined"
        << endl;
    cout << ListOfChannels()
        << endl;
    exit(0);
}

string InclusiveProcess::ListOfChannels() {
    stringstream stream;
    for (int i = 0; i < _channels.size(); i++) {
        stream << _channels[i]->Name() << endl;
    }
    return stream.str();
}




AsSeries InclusiveProcess::QCDCorrections() {
    AsSeries res(2, 0.0);
    for (int j = 0; j < _channels.size(); j++) {
        res = res+_channels[j]->QCDResult();
    }
    return res;
}

AsSeries InclusiveProcess::MurEvolution(const AsSeries& expansion,
                                        const double& L) {
    const double b0 = consts::beta_zero;
    const double b1 = consts::beta_one;
    const double b2 = consts::beta_two;
    const ResultPair c2 = expansion.term_of_order(2);
    const ResultPair c3 = expansion.term_of_order(3);
    const ResultPair c4 = expansion.term_of_order(4);
    const ResultPair c5 = expansion.term_of_order(5);
    AsSeries res2(2, c2);
    AsSeries res3(3, c3 - 2* b0* c2* L);
    AsSeries res4(4, c4 + (-2*b1*c2-3*b0*c3)*L + 3*c2*pow(b0*L, 2));
    AsSeries res5(5, c5 + (-2 * b2 * c2 - 3 * b1 * c3 - 4 * b0 * c4) * L
                  + (7 * b0 * b1*  c2 + 6 * pow(b0, 2.) *c3)*pow(L, 2.)
                  + (-4 *pow(b0, 3.)* c2)*pow(L, 3.));
    return res2+res3+res4+res5;
}

AsSeries InclusiveProcess::GenericMurEvolution(const AsSeries& expansion,
                                        const double& L) {
    const double b0 = consts::beta_zero;
    const double b1 = consts::beta_one;
    const double b2 = consts::beta_two;
    const ResultPair c0 = expansion.term_of_order(0);
    const ResultPair c1 = expansion.term_of_order(1);
    const ResultPair c2 = expansion.term_of_order(2);
    const ResultPair c3 = expansion.term_of_order(3);
    const ResultPair c4 = expansion.term_of_order(4);
    const ResultPair c5 = expansion.term_of_order(5);
    AsSeries res0(0, c0);
    AsSeries res1(1, c1);
    AsSeries res2(2, c2 - b0*c1*L);
    AsSeries res3(3, c3 - 2*b0*c2*L + c1*(-(b1*L) + pow(b0*L, 2)));
    AsSeries res4(4, c4 - b2*c1*L - 2*b1*c2*L - 3*b0*c3*L
                  + (5./2. *b0*b1*c1*pow(L, 2))
                  + 3*pow(b0, 2)*c2*pow(L, 2) - pow(b0, 3)*c1*pow(L, 3));
    AsSeries res5(5, c5 - 2*b2*c2*L - 3*b1*c3*L - 4*b0*c4*L
                  + 7*b0*b1*c2*pow(L, 2) + 6*pow(b0, 2)*c3*pow(L, 2)
                  - 4*pow(b0, 3)*c2*pow(L, 3));
    return res0+res1+res2+res3+res4+res5;
}










