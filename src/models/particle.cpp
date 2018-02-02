
#include "particle.h"
#include "iostream"
#include <stdlib.h> //: for exit()

using namespace std;



void CouplingConstant::evolve(const double& mur,int porder,const double zmass)
{
    const int N=100000;
    
    double as_prev = v_ / consts::Pi;
    double step=log(pow(mur,2.0)/pow(zmass,2.0))/N;
    if (porder==0)
        {
        for (int i=0;i<N;i++)
            {
            double incr = step * pow(as_prev,2.0) * consts::beta_zero;
            as_prev = as_prev - incr;
            }
        }
    else if (porder==1)
        {
        for (int i=0;i<N;i++)
            {
            double incr = step * pow(as_prev,2.0)
                                *(  consts::beta_zero
                                  + consts::beta_one * as_prev);
            as_prev = as_prev - incr;
            }
        }
    else if (porder==2)
        {
        for (int i=0;i<N;i++)
            {
            double incr = step * pow(as_prev,2.0)
                            *( consts::beta_zero
                             + consts::beta_one * as_prev
                             + consts::beta_two * pow(as_prev,2.0)
                              );
            as_prev = as_prev - incr;
            }
        }
    else if (porder==3)
    {
        for (int i=0;i<N;i++)
        {
            double incr = step * pow(as_prev,2.0)
                                *( consts::beta_zero
                                 + consts::beta_one * as_prev
                                 + consts::beta_two * pow(as_prev,2.0)
                                 + consts::beta_three * pow(as_prev,3.0)
                                                   );
            as_prev = as_prev - incr;
        }
    }
    else
        {
        cout<<"\n["<<__func__<<"]\t  : Unknown perturbative order : "<<porder<<endl;
        exit(0);
        }
    v_ = as_prev*consts::Pi;
}



//------------------------------------------------------------------------------

void MassParameter::evolve(CouplingConstant as,const double & mur,
                           int porder,const double& zmass)
{
    string verbosity_level="normal";//"normal";
    if (scheme_=="msbar")
        {
        // evolving as to the reference scale of this mass
        as.evolve(ref_scale_,porder,zmass);
        double as_at_ref_scale = as.v();
        //
        const int N=100000;
        
    
    
        double as_prev= as_at_ref_scale/consts::Pi;
        double m_prev = m_at_ref_scale_;
    
        double step=log(pow(mur,2.0)/pow(ref_scale_,2.0))/N;
        for (int i=0;i<N;i++)
            {
            double incr;
            double m_incr;
            switch(porder)
                {
                case 0:
                incr = step * pow(as_prev,2.0) * consts::beta_zero;
                m_incr = step * consts::gamma_one * as_prev * m_prev;
                break;
                case 1:
                incr = step * pow(as_prev,2.0)*(consts::beta_zero
                                                +consts::beta_one*as_prev);
                m_incr = step * ( consts::gamma_one*as_prev+
                                 + consts::gamma_two * pow(as_prev,2.0)
                                 )* m_prev;
                break;
                case 2:
                incr = step * pow(as_prev,2.0)*(consts::beta_zero
                                                +consts::beta_one*as_prev
                                                + consts::beta_two * pow(as_prev,2.0));
                m_incr = step * ( consts::gamma_one*as_prev
                                 + consts::gamma_two * pow(as_prev,2.0)
                                 + consts::gamma_three * pow(as_prev,3.0)
                                 )* m_prev;
                break;
                case 3:
                incr = step * pow(as_prev,2.0)*(consts::beta_zero
                                                +consts::beta_one*as_prev
                                                + consts::beta_two * pow(as_prev,2.0)
                                                + consts::beta_three * pow(as_prev,3.0));
                m_incr = step * ( consts::gamma_one*as_prev
                                + consts::gamma_two * pow(as_prev,2.0)
                                + consts::gamma_three * pow(as_prev,3.0)
                                + consts::gamma_four * pow(as_prev,4.0)
                                     )* m_prev;
                break;
                default:cout<<"\n["<<__func__<<"]\t  : Unknown perturbative order : "<<porder<<endl;
                exit(0);
                }
            as_prev = as_prev - incr;
            m_prev = m_prev - m_incr;
            }
        m_=m_prev;
            if (verbosity_level=="debugging"){
                cout << " msbar scheme : evolved from m("<<ref_scale_<<")="
                     << m_at_ref_scale_
                     << " -> m("<<mur<<")="<<m_<<endl;
            }
        }
    else
        {
        }
}

void MassParameter::evolve_to_4_loops(CouplingConstant as,
                                      const double & mur,
                                      const double& zmass)
{
    string verbosity_level="normal"; //"normal";
    if (scheme_=="msbar")
    {
        // evolving as to the reference scale of this mass
        as.evolve(ref_scale_,3,zmass);
        double as_at_ref_scale = as.v();
        //
        const int N=100000;
        double as_prev= as_at_ref_scale/consts::Pi;
        double m_prev = m_at_ref_scale_;
        
        double step=log(pow(mur,2.0)/pow(ref_scale_,2.0))/N;
        for (int i=0;i<N;i++)
        {
            double incr;
            double m_incr;
            incr = step * pow(as_prev,2.0)*(  consts::beta_zero
                                            + consts::beta_one * as_prev
                                            + consts::beta_two * pow(as_prev,2.0)
                                            + consts::beta_three * pow(as_prev,3.0)
                                            );
            m_incr = step * (  consts::gamma_one * as_prev
                             + consts::gamma_two * pow(as_prev,2.0)
                             + consts::gamma_three * pow(as_prev,3.0)
                             + consts::gamma_four * pow(as_prev,4.0)
                             )* m_prev;
            
            as_prev = as_prev - incr;
            m_prev = m_prev - m_incr;
        }
        m_=m_prev;
        if (verbosity_level=="debugging"){
            cout<<" msbar scheme : evolved from m("<<ref_scale_<<")="
            <<m_at_ref_scale_<<" -> m("<<mur<<")="<<m_<<endl;
        }
    }
    else
    {
        if (verbosity_level=="debugging"){
            cout<<" on-shell scheme, no evolution. Mass used = "<<m_<<endl;
        }
    }
}

//------------------------------------------------------------------------------
Particle::Particle() {
    name_="noname";
    Y_=0.0;
    width_=0.0;
    charge_ = 0.0;
    m_ = NULL;
}


Particle::Particle(const string& name)
{
    name_=name;
    Y_=0.0;
    width_=0.0;
    charge_ = 0.0;
    m_ = NULL;
}

void Particle::set_msbar_mass(const double& m,const double& ref_scale)
{
    m_ = new MassParameter(m,ref_scale);
}

void Particle::set_pole_mass(const double & m)
{
    m_ = new MassParameter(m);
}

void Particle::consolidate(CouplingConstant as,const double& mur, int porder,
                                 const double & mh,const double& zmass)
{
   if (m_ == NULL)
       {
       cout<<"Error in Model/Particle : mass has not been initialized before trying to evolve it"<<endl<<endl;
       exit(0);
       }
    else
        {
            
            m_->evolve(as,mur,porder,zmass);
            _m_higgs = mh;
        }
}





complex<double> Particle::cm_sq(){
    return complex<double>(pow(m_->value(),2.0), -width_*m_->value());
}

complex<double> Particle::X(){return - Wq()/pow(sqrt(1.0-Wq())+1.0,2.0);}
complex<double> Particle::Wq(){return 4.0*cm_sq()/pow(_m_higgs,2.);}
















