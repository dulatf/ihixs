#include "gluon_fusion_ew_coefficients.h"

#include "nlo_exact_matrix_elements.h"


//------------------------------------------------------------------------------

GluonFusionEWCoefficients::GluonFusionEWCoefficients(const CModel& model)
{
    
#include "electroweak_data.h"
    
    // reads EWK  corrections from Fig. 21 of
    //      http://arXiv.org/pdf/0809.3667
    //  by Actis, Passarino, Sturm, Uccirati
    //   Data in  file "./electroweak.h" provided by the authors.
    //   They have chosen the following paramegters:
    
    //       Mw = 80.398
    //       GammaW = 2.093
    //       Mz = 91.1876
    //       GammaZ = 2.4952
    //       Gfermi = 1.16637e-5 !/GeV^2
    //       a(0) = 1.0/137.0359911.0
    //       alphas_Mz = 0.118.0
    //       Mtop =  170.9
    
    // note that we are *not* using the msbar top mass
    // that the rest of the program might use, here, to make sure that the eff_mh
    // is indeed mh * 170.9/current_pole_mass_of_top
    // This induces a per mille difference with ihixs (for low mh)
    // ehixs = ihixs - 0.1%
    //model_=model;
    
    //cout<<"\nCalcualting ew correction factor";
    const double mtop_pass = 170.9;
    const double current_mt_mass_os = 172.5;
	const double eff_mh =model.higgs.m()*mtop_pass/current_mt_mass_os;
    
    if (eff_mh<100.0 or eff_mh>500.0)
    {
        cout<<"\nSoft ew corrections not available (mh<100 or mh>500). They are set to zero since mh = "<<eff_mh
        <<" mtop_pass/ mtop = "<< mtop_pass/current_mt_mass_os<<endl;
        cout<<" mh from model = "<<model.higgs.m();
        NLO_ew_coeff_ = 0.0;
    }
    else
    {
        int N = ew_data.size();
        int position=0;
        
        for (int i=0;i<ew_data.size()-1;i++)
        {
            const double mleft = ew_data[i]->mass;
            const double mright = ew_data[i+1]->mass;
            if (mleft < eff_mh and eff_mh<mright)
            {
                position = i;
                if (mright-eff_mh < eff_mh-mleft) {position = i+1;}
                break;
            }
        }
        int ibefore,ihere,iafter;
        if (position==0){ibefore = 0;ihere=1;iafter=2;}
        else if (position==N-1){ibefore = N-3;ihere=N-2;iafter=N-1;}
        else
        {
            
            ibefore = position-1;
            ihere = position;
            iafter = position+1;
        }
        double x[3]={ew_data[ibefore]->mass,
            ew_data[ihere]->mass,
            ew_data[iafter]->mass};
        double y[3]={ew_data[ibefore]->deltaew,
            ew_data[ihere]->deltaew,
            ew_data[iafter]->deltaew};
        vector<double> c = givecoeff(x,y);
        double res = (c[2] + c[1] * eff_mh + c[0] * eff_mh*eff_mh)/100.0;
        // calculating born for top only with mtop = 170.9 GeV
        // which was the choice of the authors of 0809.3667
        complex<double> Wt(4.0*170.9*170.9/pow(eff_mh,2.0),0.0);
        complex<double> xt = -Wt/pow(sqrt(1.0-Wt)+1.0,2.0);
        complex<double> born_special = h_exact::born(xt);
        NLO_ew_coeff_ = (sqrt(1.0+res)-1.0 ) * born_special;
    }
}


vector<double>  GluonFusionEWCoefficients::givecoeff(double x[3],double y[3])
{
    vector<double> res;
    const double dx12 = x[0]-x[1];
    const double dx23 = x[1]-x[2];
    const double dx31 = x[2]-x[0];
    
    const double den=dx12*dx23*dx31;
    res.push_back((-y[0]*dx23-y[1]*dx31-y[2]*dx12)/den);
    
    res.push_back(( y[0]*(x[1]*x[1]-x[2]*x[2])
                   +y[1]*(x[2]*x[2]-x[0]*x[0])
                   +y[2]*(x[0]*x[0]-x[1]*x[1]) ) /den);
    
    res.push_back( (-y[0]*dx23*x[1]*x[2]
                    -y[1]*dx31*x[0]*x[2]
                    -y[2]*dx12*x[0]*x[1]) / den);
    return res;
}






//------------------------------------------------------------------------------


