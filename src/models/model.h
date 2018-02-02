#ifndef MODEL_H
#define MODEL_H

#include "particle.h"
#include<vector>
#include <stdlib.h> //: for exit()
#include "user_interface.h"//: for UI
using namespace std;



class CModel{
public:
     CModel();
     ~CModel(){};
         
    void Configure(const double& a_at_mz,const double & mur_over_mh,
                      int porder,const double& mh);
    double alpha_strong(){return alpha_s->v();}

    double mu_r() const {return mu_r_;}
    
    Particle higgs;
    vector<Particle*> quarks;
    vector<Particle*> leptons;
    vector<Particle*> vector_bosons;
    Particle top;
    VectorBoson W;
    VectorBoson Z;
    
    Particle bottom;
    Particle charm;
    
    void ReadParameters(const UserInterface& UI);
    
    
private:
	double mu_r_;
    double _porder;
    CouplingConstant* alpha_s;
    
    CModel(const CModel&);// Prevent copy construction
    CModel& operator=(const CModel&);// Prevent assignment
    
private:
    void check_that_particle_masses_are_set();
    

    
     

};


#endif
