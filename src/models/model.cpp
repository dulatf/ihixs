

#include "model.h"
#include "iostream"

using namespace::std;


CModel::CModel()
{
    // quark masses set by user!
    
    top=Particle("top");
    top.set_charge(2.0/3.0);
    top.set_Y(1.0);
    top.set_width(0.0);
    //
    bottom=Particle("bottom");
    bottom.set_charge(-1.0/3.0);
    bottom.set_Y(1.0);
    bottom.set_width(0.0);
    //
    charm=Particle("charm");
    charm.set_charge(2.0/3.0);
    charm.set_Y(1.0);
    charm.set_width(0.0);
    //
    W=VectorBoson("W");
    W.set_pole_mass(80.403);
    W.set_charge(1.0);
    W.set_Y(1.0);
    W.set_width(2.141);
    
    //
    Z=VectorBoson("Z");
    Z.set_pole_mass(91.1876);
    Z.set_charge(0.0);
    Z.set_Y(1.0);
    Z.set_width(2.4952);
    //
    //
    higgs=Particle("higgs");
    higgs.set_pole_mass(125.0);
    higgs.set_charge(0.0);
    higgs.set_Y(1.0);
    higgs.set_width(4.36E-003);
    //
    quarks.push_back(&top);
    quarks.push_back(&bottom);
    quarks.push_back(&charm);
    vector_bosons.push_back(&W);
    vector_bosons.push_back(&Z);

    mu_r_=0.0;
}


void CModel::ReadParameters(const UserInterface& UI){
    // quark masses
    if (UI.giveString("top_scheme")=="msbar")
        top.set_msbar_mass(UI.giveDouble("mt_msbar"),
                           UI.giveDouble("mt_msbar_ref_scale"));
    else if (UI.giveString("top_scheme")=="on-shell")
        top.set_pole_mass(UI.giveDouble("mt_on_shell"));
    else{
        cout<<"Model: unrecognized scheme for top "<<UI.giveString("top_scheme")<<endl;
        exit(EXIT_FAILURE);
    }
    
    if (UI.giveString("bottom_scheme")=="msbar")
        bottom.set_msbar_mass(UI.giveDouble("mb_msbar"),
                            UI.giveDouble("mb_msbar_ref_scale"));
    else if (UI.giveString("bottom_scheme")=="on-shell")
                    bottom.set_pole_mass(UI.giveDouble("mb_on_shell"));
    else{
        cout<<"Model: unrecognized scheme for bottom "<<UI.giveString("bottom_scheme")<<endl;
        exit(EXIT_FAILURE);
    }
    
    if (UI.giveString("charm_scheme")=="msbar")
        charm.set_msbar_mass(UI.giveDouble("mc_msbar"),UI.giveDouble("mc_msbar_ref_scale"));
    else if (UI.giveString("charm_scheme")=="on-shell") charm.set_pole_mass(UI.giveDouble("mc_on_shell"));
    else{
        cout<<"Model: unrecognized scheme for charm "<<UI.giveString("charm_scheme")<<endl;
        exit(EXIT_FAILURE);
    }
    
    top.set_Y(UI.giveDouble("y_top"));
    bottom.set_Y(UI.giveDouble("y_bot"));
    charm.set_Y(UI.giveDouble("y_charm"));
    
    top.set_width(UI.giveDouble("gamma_top"));
    bottom.set_width(UI.giveDouble("gamma_bot"));
    charm.set_width(UI.giveDouble("gamma_charm"));
}



void CModel::check_that_particle_masses_are_set(){
    for (int i=0;i<quarks.size();i++){
        if (quarks[i]->MassIsNotSet()){
            cout<<"[model]: Error: the mass of "<<quarks[i]->name()
            <<" was not set properly "<<endl<<endl;
            exit(EXIT_FAILURE);
        }
    }
}


//: remove the ratio mur_over_mh in Configure
void CModel::Configure(const double & a_at_mz,
                         const double & mur_over_mh, int porder,
                         const double& mh)
{
    
    string verbosity_level="silent";//"silent";//"debugging";
    check_that_particle_masses_are_set();
    higgs.set_pole_mass(mh);
    
    const double new_mur = mur_over_mh * higgs.m();
    if ((new_mur != mu_r_) or (porder != _porder))
    {
        mu_r_ = mur_over_mh * higgs.m();
        if (porder<0 or porder>3){
            cout<<"[model] Error: porder fed to model is not 0,1,2 or 3: porder="<<porder<<endl;
            exit(EXIT_FAILURE);
        }
        
        
        else _porder = porder;

        if (verbosity_level=="debugging"){
            cout<<endl<<"[model] couplings and masses evolved to order "<<_porder<<endl;
            }
        
        alpha_s = new CouplingConstant(a_at_mz,Z.m());
        
        for (int i=0;i<quarks.size();i++) {
            quarks[i]->consolidate(*alpha_s,mu_r_,_porder,higgs.m(),Z.m());
        }
        for (int i=0;i<vector_bosons.size();i++)
            vector_bosons[i]->consolidate(*alpha_s,mu_r_,_porder,higgs.m(),Z.m());
        alpha_s->evolve(mu_r_,_porder,Z.m());
        if (verbosity_level=="debugging"){
            cout<<"[model] a_s(m_z) = "<<a_at_mz<<" -> a_s("<<mu_r_<<") = "<<
                alpha_s->v()<<endl;
            }
        
        // complex masses in gw,sw,cw
        complex<double> gw = sqrt(4.0*sqrt(2.0)*consts::G_fermi * W.cm_sq());
        complex<double> cw_sq = W.cm_sq()/Z.cm_sq();
        complex<double> cw = sqrt(cw_sq);
        complex<double> sw_sq = 1.0-cw_sq;
        //real masses in gw, sw, cw
        gw = sqrt(4.0*sqrt(2.0)*consts::G_fermi * pow(W.m(),2.0));
        cw_sq = pow(W.m()/Z.m(),2.0);
        cw = sqrt(cw_sq);
        sw_sq = 1.0-cw_sq;
        
        W.cv_up = gw / sqrt(2.0);
        W.cv_down = gw / sqrt(2.0);
        W.ca_up = gw / sqrt(2.0);
        W.ca_down = gw / sqrt(2.0);
        W.lamda = 2.0;
        
        
        Z.cv_up = gw/cw*(1.0/2.0 - 4.0/3.0 * sw_sq);
        Z.cv_down = gw/cw*(-1.0/2.0 + 2.0/3.0 * sw_sq);
        Z.ca_up = gw/cw*1.0/2.0;
        Z.ca_down = gw/cw*1.0/2.0;
        Z.lamda = 2.0;
        if (verbosity_level=="debugging"){
            cout<<"[model] gw = "<<gw<<"\tsw^2 = "<< sw_sq
                <<"\tcw = "<< cw<<endl;
            cout<<"[model] G_F = "<<consts::G_fermi<<endl;
            cout<<"[model] mw = "<<W.m()<<endl;

            cout<<"[model] gw_up for W = "<<pow(W.cv_up,2.0) + pow(W.ca_up,2.0)<<endl;
            cout<<"[model] gw_up for Z = "<<pow(Z.cv_up,2.0) + pow(Z.ca_up,2.0)<<endl;
            cout<<"[model] lambda for W = "<<W.lamda<<endl;
            cout<<"[model] lambda for Z = "<<Z.lamda<<endl;
        }
    }
    
}



























