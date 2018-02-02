// Copyright 2015 ihixs team

#ifndef SRC_CORE_INCLUSIVE_PROCESS_H_
#define SRC_CORE_INCLUSIVE_PROCESS_H_

#include <vector>
#include <string>
#include "src/tools/user_interface.h"
#include "src/core/channel.h"
#include "src/models/model.h"

class InclusiveProcess{
 public:
    // constructor
    explicit InclusiveProcess(const UserInterface& UI);
    void Initialize(const UserInterface& UI);
    ~InclusiveProcess();
    // utility function returning the name of the process
    string Name() const {return _name;}
    // contributions determined by daughters
    virtual void SetUpContributions() = 0;
    // post vegas evaluation depending on particular process
    virtual void EvaluatePostVegasDependingOnProcess() = 0;
    // evaluates all terms defined
    void Evaluate();
    // evaluates single term
    void EvaluateSingleTerm(const string& term_name);
    // prints out all terms defined
    friend ostream& operator<<(ostream&, const InclusiveProcess&);
    // prints out all channel names
    string ListOfChannels();
    // channel breakdown (of whatever channel is defined)
    // REFACTORING: this is probably useless
    string ChannelBreakdown();
    // QCD corrections as a series in as
    // REFACTORING: what about ew corrections?
    // should this be in daughter classes?
    AsSeries QCDCorrections();
    InputParameters *GetInput() { return &_input;}
    //
    string InputInformation() {return _input.InputInformation();}
    //
    Channel* GiveChannel(const string& name);
    //
    double AsOverPi() {return _input._as_over_pi;}

 protected:
    // all input for the xs
    InputParameters _input;
    // all terms defined
    vector<Channel*> _channels;
    string _name;

 protected:
    AsSeries MurEvolution(const AsSeries& expansion, const double& L);
    AsSeries GenericMurEvolution(const AsSeries& expansion, const double& L);

 private:
    void Truncate();
    int DetermineQCDPerturbativeOrder(const string&);
    vector<string> DetermineChannels();
    bool ChannelNameExists(const string& name,
                           const vector<string>& channel_list);
};




#endif  // SRC_CORE_INCLUSIVE_PROCESS_H_
