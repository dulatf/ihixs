#ifndef PARTICLEOBJECT_H
#define PARTICLEOBJECT_H

#include<complex>
#include<string>
#include "constants.h"
#include <iostream>
using namespace std;



class CouplingConstant{
public://methods
    CouplingConstant(double v,const double& ref_scale): v_(v),_name("noname"){};
    CouplingConstant(double v,const string& name)
    : v_(v),_name(name){};
    double v(){return v_;}
    void evolve(const double& mu,int porder,const double zmass);
    string name(){return _name;}
private://data
    double v_;
    string _name;
    //double ref_scale_;
};


class MassParameter{
public://methods
    MassParameter(double m,const double ref_scale){scheme_ = "msbar"; m_=0; ref_scale_ = ref_scale; m_at_ref_scale_=m;}
    MassParameter(double m){m_ = m;m_at_ref_scale_=m; ref_scale_ = 0.0; scheme_ = "on-shell";}
    double value(){return m_;}
    void evolve(CouplingConstant as,const double& mu,int porder,const double& zmass);
    void evolve_to_4_loops(CouplingConstant as,
                                          const double & mur,
                                          const double& zmass);
    string scheme(){return scheme_;}
    
    bool MassIsNotSet(){
        if (
                m_at_ref_scale_<1e-5
                or
                (scheme_ == "msbar" and ref_scale_<0.1)
            )
            return true;
        else
            return false;
            }
private://data
    double m_;
    double ref_scale_;
    double m_at_ref_scale_;
    string scheme_;
};

class Particle
{
public:
    Particle();
    ~Particle() {};
    Particle(const string&);//: name
    //
    void set_msbar_mass(const double& mass,const double & scale);
    void set_pole_mass(const double& mass);
    void set_width(const double& w){width_=w;}
    void set_charge(const double& ch){charge_ = ch;}
    void set_Y(const double& Y){Y_ = Y;}
    //
    void consolidate(CouplingConstant as,const double& mur, int porder,
                     const double& mh,const double& zmass);
    //
    double width(){return width_;}
    double m()const {return m_->value();}
    complex<double> cm_sq();
    double Y(){return Y_;}
    complex<double> X();
    complex<double> Wq();
    double charge(){return charge_;}
    string name(){return name_;}
    string scheme(){return m_->scheme();}
    
    bool MassIsNotSet(){return m_->MassIsNotSet();}
    
private://data
    string name_;
    double width_;
    double charge_;
    double Y_;
    MassParameter* m_;
	
    double _m_higgs;
    double m_at_scale_mur;
    double evolved_scale;
    bool evolved;
    bool needs_evolution;
private://methods
    double evolution_step();


};

class VectorBoson: public Particle
{
public:
    VectorBoson(){};
    VectorBoson(const string& name): Particle(name){};
    complex<double> cv_up;
    complex<double> cv_down;
    complex<double> ca_up;
    complex<double> ca_down;
    complex<double> lamda;
private:

};



#endif









