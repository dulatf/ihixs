

    
#include<complex>
#include "constants.h"
#include "cppchaplin.h"

#include "model.h"
#include "nlo_exact_matrix_elements.h"
#include <iostream>
#include <iomanip>
using namespace std;



namespace h_exact {
#define complex_double std::complex<double>


    double ExactRegularCalculator::gg_reg_Log0(
                    const double& z,
                    const double& lambda) {
        
        return    gg_reg_Log0_hard(z,lambda)
                + gg_reg_Log0_collinear(z,lambda);
    }
    
    double ExactRegularCalculator::gg_reg_Log0_hard(
                                               const double& z,
                                               const double& lambda) {
        
        const double H = sum_of_abs_sq_of_Aqi(z,lambda);
        
        const double limit = pow(1.- z + z*z,2.) * _B_sq;
        
        const double hard = 3.* (1./2. * pow(z,4.) * H - limit)
        /z / (1.-z) / lambda / (1.-lambda);
        
        return hard;
    }
    
    double ExactRegularCalculator::gg_reg_Log0_collinear(
                                               const double& z,
                                               const double& lambda) {
        
        
        const double pgg = 1./z + z*(1.-z)-2.;
        
        const double collinear = _B_sq * (6. * pgg * log(pow(1.-z,2)/z)
                                          -6.*log(z)/(1.-z)
                                          );
        
        
        
        return collinear;
    }

    double ExactRegularCalculator::gg_reg_Log1(
                                               const double& z,
                                               const double& lambda)
    {
        
        const double pgg = 1./z + z*(1.-z)-2.;
                                
        const double collinear = _B_sq * (- 6. * pgg);
        return collinear;
    }
    
    double ExactRegularCalculator::qg_reg_Log0(const double& z,const double& lambda)
    {
        return  qg_reg_Log0_hard(z,lambda)
                +qg_reg_Log0_collinear(z,lambda);
    }
    
    double ExactRegularCalculator::qg_reg_Log0_hard(const double& z,const double& lambda)
    {
        const double ysp = -(1.-z)*(1.-lambda)/z;
        const double sum_Aqqgh_lambda=sum_of_abs_sq_of_Aqqgh(1./ysp);
        
        const double term2 = (
                              sum_Aqqgh_lambda * (1.+pow(1.-z,2.)*pow(lambda,2.))/z
                              - _B_sq * (1.+pow(1.-z,2.))/z
                              )
        /(1.-lambda);
        
        return (2./3.*term2);
    }
    
    double ExactRegularCalculator::qg_reg_Log0_collinear(const double& z,const double& lambda)
    {
        const double Cf= 4./3.;
        const double Pgq = Cf/2.0*(1.0+pow(1.0-z,2.))/z;
        const double term1 = _B_sq * (Cf/2.*z
                                      - Pgq * log(z/pow(1.-z,2.))
                                      );
        return (term1);
    }
    
    double ExactRegularCalculator::qg_reg_Log1(const double& z,const double& lambda)
    {
        const double Cf= 4./3.;
        const double Pgq = Cf/2.0*(1.0+pow(1.0-z,2.))/z;
        const double term1 = _B_sq * (- Pgq ); //L=log(muf**2/mH**2)
        return (term1);
    }
    
    double ExactRegularCalculator::qqb_reg_Log0(const double& z,const double& lambda)
    {
        return  32./27. * pow(1.-z, 3.) /z * sum_of_abs_sq_of_Aqqgh(z);
    }
    
    
    complex<double> ExactRegularCalculator::sum_of_Aqqgh(const double& y)
    {
        complex<double> res(0.0,0.0);
        for (int i=0;i<_factor.size();i++){
           res += h_exact::Aqqgh_cpp(y,_Xfactor[i]) * _factorqqgh[i];
        }
        return res;
    }
    
    
    double ExactRegularCalculator::sum_of_abs_sq_of_Aqqgh(const double& y)
    {
        return pow(abs(sum_of_Aqqgh(y)),2.0);
    }
    
    

    
    void ExactRegularCalculator::compute_soft()
    {
        
        complex<double> V(0.0,0.0) ;
        complex<double> Born(0.0,0.0);
        string verbosity_level="non-debugging";
        
        for (int i=0;i<_model->quarks.size();i++)
        {
            // every quark has a scheme: on-shell or msbar
            // the ggf_exact_virtual_ep0 below is equal to
            // A + B * scheme_dependent_coeff
            // where scheme_dependent_coeff = 4/3 for on-shell scheme
            // and   scheme_dependent_coeff = log(mq^2/mur^2) for MS_bar
            complex<double> scheme_dependent_coeff;
            if (_model->quarks[i]->scheme()=="on-shell")
                scheme_dependent_coeff = complex<double> (4.0/3.0,0.0);
            else if (_model->quarks[i]->scheme()=="msbar")
            {
                scheme_dependent_coeff =  log(_model->quarks[i]->cm_sq()
                                              / pow(_model->mu_r(),2.)
                                              );
            }
            V = V + _model->quarks[i]->Y() * h_exact::ggf_exact_virtual_ep0(_model->quarks[i]->X(),scheme_dependent_coeff);
            Born = Born+_model->quarks[i]->Y() * born(_model->quarks[i]->X());
            if (verbosity_level=="debugging"){
                cout<<"[GluonFusionExactCoefficients] NLO "<<_model->quarks[i]->name()
                <<setprecision(4)
                <<" : contributes "<<_model->quarks[i]->Y() * h_exact::ggf_exact_virtual_ep0(_model->quarks[i]->X(),scheme_dependent_coeff)
                <<"\t x= "<<_model->quarks[i]->X()<<"\t W = "<<_model->quarks[i]->Wq()
                <<"\t m = "<<_model->quarks[i]->m()
                <<"\t cm_sq = "<<_model->quarks[i]->cm_sq()
                <<endl;
            }
        }
        
        double res = 2.0 * real(Born * conj(V))
        + pow(abs(Born),2.0) * consts::pi_square;
        if (verbosity_level=="debugging"){
            cout<<"\n NLO soft= "<<res<<"\t"<<conj(V);
            cout<<"\n NLO soft : 2.0*real(B*Vstar) = "<<2.0*real(Born * conj(V))
            <<" |B|^2*pi^2 = "<<pow(abs(Born),2.0) * consts::pi_square;
        }
        _soft = res;
    }
    
  void ExactRegularCalculator:: SetUp() {
      
      _mh = _model->higgs.m();
      for (int i=0;i<_model->quarks.size();i++){
          double Yq = _model->quarks[i]->Y();
         
          if (Yq!=0.0)
          {
              
              complex<double> Wq = _model->quarks[i]->Wq();
              complex<double> mq_cplx = sqrt(
                                             _model->quarks[i]->cm_sq());
              complex<double> factor = Yq * Wq * 3.0/32.0 /pow(_mh,2.0);
              complex<double> factorqqgh = Yq * Wq * 3.0 / 4.0;

              _factor.push_back(factor);
              _factorqqgh.push_back(factorqqgh);
              _mq_cplx.push_back(mq_cplx);
              _Xfactor.push_back(_model->quarks[i]->X());
          }
      }
      _B_sq = pow(abs(h_exact::born_exact_summed_over_quarks(*_model)),2.);

  }
  
    double ExactRegularCalculator::sum_of_abs_sq_of_Aqi(const double &z,const double & lambda)
    {
        complex<double> A1(0.0,0.0);
        complex<double> A2(0.0,0.0);
        complex<double> A3(0.0,0.0);
        complex<double> A4(0.0,0.0);
        for (int i=0;i<_factor.size();i++) {
                A1 += h_exact::Aq1(z,lambda,_mq_cplx[i],_mh) * _factor[i];
                A2 += h_exact::Aq2a(z,lambda,_mq_cplx[i],_mh) * _factor[i];
                A3 += h_exact::Aq2b(z,lambda,_mq_cplx[i],_mh) * _factor[i];
                A4 += h_exact::Aq2c(z,lambda,_mq_cplx[i],_mh) * _factor[i];
        }
        return(    pow(abs(A1),2.0)
               + pow(abs(A2),2.0)
               + pow(abs(A3),2.0)
               + pow(abs(A4),2.0));
    }
    

    
    
    
    int point_counter=0;
    
    double gg_reg_exact_nlo_L0(const double& z,
                               const double& lambda,
                               const CModel& model)
    {
        //L=log(muf^2 / mt^2)
        //     regular part of NLO Real corrections with massive quarks
        //
        //      ggreg_nlo = 1/2 * Nc * as/Pi *
        //               [ z^4  S_i |S_q real_ggg_ampl_i(m_q)|^2
        //              - (1-z+z^2) * 2 * |S_q born_q|^2 ]
        //              --------------------------------------------------------
        //                 z * (1-z) * lambda * (1-lambda)
        
        //    sum_over_quarks: sums over all massive quarks the function
        //    that is passed to it as its first argument and it also
        //    multiplies by a normalization factor Wq(i)*3.0/32d0/QQQ**2
        //    which is tau_q*3/32/m_higgs^2
        
        const double H = sum_of_abs_sq_of_Aqi(z,lambda,model);
        const double B_sq = pow(abs(born_exact_summed_over_quarks(model)),2.);
        const double limit = pow(1.- z + z*z,2.) * B_sq;
        const double pgg = 1./z + z*(1.-z)-2.;
        
        const double hard = 3.* (1./2. * pow(z,4.) * H - limit)
        /z / (1.-z) / lambda / (1.-lambda);
        
        const double collinear = B_sq * (6. * pgg * log(pow(1.-z,2)/z)
                                         -6.*log(z)/(1.-z)
                                         );
        point_counter++;
        return hard + collinear;
        
    }
    
    double gg_reg_exact_nlo_L0_hard_only(const double& z,
                               const double& lambda,
                               const CModel& model)
    {
        //L=log(muf^2 / mt^2)
        //     regular part of NLO Real corrections with massive quarks
        //
        //      ggreg_nlo = 1/2 * Nc * as/Pi *
        //               [ z^4  S_i |S_q real_ggg_ampl_i(m_q)|^2
        //              - (1-z+z^2) * 2 * |S_q born_q|^2 ]
        //              --------------------------------------------------------
        //                 z * (1-z) * lambda * (1-lambda)
        
        //    sum_over_quarks: sums over all massive quarks the function
        //    that is passed to it as its first argument and it also
        //    multiplies by a normalization factor Wq(i)*3.0/32d0/QQQ**2
        //    which is tau_q*3/32/m_higgs^2
        
        const double H = sum_of_abs_sq_of_Aqi(z,lambda,model);
        const double B_sq = pow(abs(born_exact_summed_over_quarks(model)),2.);
        const double limit = pow(1.- z + z*z,2.) * B_sq;
        //const double pgg = 1./z + z*(1.-z)-2.;
        
        const double hard = 3.* (1./2. * pow(z,4.) * H - limit)
        /z / (1.-z) / lambda / (1.-lambda);
        
        // we remove the integrated (over lambda) collinear counterterm here
        //const double collinear = B_sq * (6. * pgg * log(pow(1.-z,2)/z)
        //                                 -6.*log(z)/(1.-z)
        //                                 );
        //point_counter++;
        //cout<<"\n"<<point_counter<<"\t"<<z<<"\t\t"<<lambda<<"\t\t"<<hard;
        return hard ;//+ collinear;
        
    }
    
    double gg_reg_exact_nlo_L1(const double& z, const double& lambda,const CModel& model)
    {
        const double B_sq = pow(abs(born_exact_summed_over_quarks(model)),2.);
        const double pgg = 1./z + z*(1.-z)-2.;
        
        const double collinear = B_sq * (- 6. * pgg);
        return collinear;
        
    }

    double gg_reg_exact_nlo(const double& z, const double& lambda,const double& L,const CModel& model)
    {
        //L=log(muf^2 / mt^2)
        //     regular part of NLO Real corrections with massive quarks
        //
        //      ggreg_nlo = 1/2 * Nc * as/Pi *
        //               [ z^4  S_i |S_q real_ggg_ampl_i(m_q)|^2
        //              - (1-z+z^2) * 2 * |S_q born_q|^2 ]
        //              --------------------------------------------------------
        //                 z * (1-z) * lambda * (1-lambda)
        
        //    sum_over_quarks: sums over all massive quarks the function
        //    that is passed to it as its first argument and it also
        //    multiplies by a normalization factor Wq(i)*3.0/32d0/QQQ**2
        //    which is tau_q*3/32/m_higgs^2
        
        const double H = sum_of_abs_sq_of_Aqi(z,lambda,model);
        const double B_sq = pow(abs(born_exact_summed_over_quarks(model)),2.);
        const double limit = pow(1.- z + z*z,2.) * B_sq;
        const double pgg = 1./z + z*(1.-z)-2.;
        
        const double hard = 3.* (1./2. * pow(z,4.) * H - limit)
                            /z / (1.-z) / lambda / (1.-lambda);
        
        const double collinear = B_sq * (6. * pgg * log(pow(1.-z,2)/z)
                                         -6.*log(z)/(1.-z)
                                         - 6.*L * pgg);
        
        return hard + collinear;
        
    }

    double qg_reg_exact_nlo(const double& z,const double& lambda,const double& L,
                            const CModel& model)
    {
        const double Cf= 4./3.;
        const double Pgq = Cf/2.0*(1.0+pow(1.0-z,2.))/z;
        const double B_sq = pow(abs(born_exact_summed_over_quarks(model)),2.);

        const double ysp = -(1.-z)*(1.-lambda)/z;
        const double sum_Aqqgh_lambda=sum_of_abs_sq_of_Aqqgh(1./ysp,model);
    
        const double term1 = B_sq * (Cf/2.*z
                                     - Pgq * log(z/pow(1.-z,2.))
                                     - Pgq * L); //L=log(muf**2/mH**2)
        const double term2 = (
             			sum_Aqqgh_lambda * (1.+pow(1.-z,2.)*pow(lambda,2.))/z
                        - B_sq * (1.+pow(1.-z,2.))/z
                              )
                            /(1.-lambda);

        return (term1+2./3.*term2);
    }
    
    
    double qg_reg_exact_nlo_L0(const double& z,const double& lambda,
                            const CModel& model)
    {
        const double Cf= 4./3.;
        const double Pgq = Cf/2.0*(1.0+pow(1.0-z,2.))/z;
        const double B_sq = pow(abs(born_exact_summed_over_quarks(model)),2.);
        
        const double ysp = -(1.-z)*(1.-lambda)/z;
        const double sum_Aqqgh_lambda=sum_of_abs_sq_of_Aqqgh(1./ysp,model);
        
        const double term1 = B_sq * (Cf/2.*z
                                     - Pgq * log(z/pow(1.-z,2.))
                                     );
        
        const double term2 = (
                              sum_Aqqgh_lambda * (1.+pow(1.-z,2.)*pow(lambda,2.))/z
                              - B_sq * (1.+pow(1.-z,2.))/z
                              )
        /(1.-lambda);
        
        return (term1+2./3.*term2);
    }
    
    
    double qg_reg_exact_nlo_L0_hard_only(const double& z,const double& lambda,
                               const CModel& model)
    {

        const double B_sq = pow(abs(born_exact_summed_over_quarks(model)),2.);
        
        const double ysp = -(1.-z)*(1.-lambda)/z;
        const double sum_Aqqgh_lambda=sum_of_abs_sq_of_Aqqgh(1./ysp,model);
        
        
        const double term2 = (
                              sum_Aqqgh_lambda * (1.+pow(1.-z,2.)*pow(lambda,2.))/z
                              - B_sq * (1.+pow(1.-z,2.))/z
                              )
        /(1.-lambda);
        
        
        // we do not add term 1 (integrated collinear counterterm) here, since the function is hard only
        return (2./3.*term2);
    }
    
    
    double qg_reg_exact_nlo_L1(const double& z,const double& lambda,
                            const CModel& model)
    {
        const double Cf= 4./3.;
        const double Pgq = Cf/2.0*(1.0+pow(1.0-z,2.))/z;
        const double B_sq = pow(abs(born_exact_summed_over_quarks(model)),2.);
        
        
        const double term1 = B_sq * (- Pgq ); //L=log(muf**2/mH**2)
        
        return (term1);
    }
    
    
    double qqb_reg_exact_nlo(const double& z,const double& lambda,const double& L,
                             const CModel& model)
    {
        return  32./27. * pow(1.-z, 3.) /z * sum_of_abs_sq_of_Aqqgh(z,model);
    }

    double qqb_reg_exact_nlo_L0(const double& z,const double& lambda,
                             const CModel& model)
    {
        return  32./27. * pow(1.-z, 3.) /z * sum_of_abs_sq_of_Aqqgh(z,model);
    }

    complex<double> F2lb(const complex<double> & x)
    {
        return 6.0*x*(pow(x,2.0)+6.0*x+1.0)*chaplin::HPL(0,0,x)
        /pow(-1.0+x,4.0)
        -6.0*x*(1.0+x)*log(x)
        /pow(-1.0+x,3.0)
        -12.0*x
        /pow(-1.0+x,2.0);
    }


    complex<double> ggf_exact_virtual_ep0(const complex<double> & x,
                                          const complex<double>& scheme_dependent_coeff)
    {
        //	Below is an implementation of the two_loop virtual
        //  amplitude in the presence of heavy quarks (the non-pi^2 piece)
        //	The implementation (by Stephan) follows eq.26-30
        //	in http://arxiv.org/abs/hep-ph/0611266
        //  see eq.28 of  http://arxiv.org/abs/hep-ph/0611266
        //
        //  It is also equal to B1*F0 of Harlander & Kant, see
        //      http://arxiv.org/abs/hep-ph/0509189
        //
        //
        //  note:   the function actually returns a real
        //          (the conjugate is added at the end)
        //  Limits:
        //	1)Its limit as m->infty is 11/2*born
        //
        //  this limit is tested by running the exact_virtual_test
        //  note:   the implementation here, from Bonciani et al (0611266)
        //          appears to have numerical instabilities as mt->infinity.
        //          In particular, the function below approaches 11/2 if mt<3TeV
        //          smoothly (it is then equal to 11/2 up to 0.0003)
        //          but the it deviates from 11/2 presumably due to cancelations
        //          between HPLs as x->1 from above on the complex circle
        //
        //  note:   the limit does not depend on the imaginary part of mq
        //          hence 11/2 is NOT a test of the correctness
        //          of the imaginary part!
        //
        //  2)mt->0 has the limit zero
        //
        //  3)mt->mh/2 is the threshold limit
        //      it is equal to 47/2-(3/2)*Pi^2+(7/4)*Pi^2*ln(2)-(49/8)*Zeta(3)
        //                      = 13.30489968
        //  all these limits are implemented as tests at exact_virtual_test.cpp
        
        complex<double> H2 = 4.0/5.0*pow(consts::z2,2.0)
                            + 2.0*consts::z3
                            + 3.0*consts::z3/2.0*log(x)
                            - 3.0*consts::z3*log(1.0-x)
                            + consts::z2*chaplin::HPL(1,0,x)
                            + 1.0/4.0*(1.0+2*consts::z2)*chaplin::HPL(0,0,x)
                            - 2.0*chaplin::HPL(1,0,0,x)
                            + chaplin::HPL(0,0,-1,0,x)
                            + 1.0/4.0*chaplin::HPL(0,0,0,0,x) + 2.0*chaplin::HPL(1,0,-1,0,x) - chaplin::HPL(1,0,0,0,x);
        
        
        complex<double> G2lCA_onehalf =
                    4.0* x/pow(x-1.0,2.0) *
                    (
                     3.0
                     + x*(1.0+8.0*x+3.0*x*x) /pow(x-1.0,3.0)*chaplin::HPL(0,0,0,x)
                     -2.0*pow(1.0+x,2.0)/pow(1.0-x,2.0)* H2
                     +consts::z3
                     -chaplin::HPL(1,0,0,x)
                     )
                    ;
        // see eq.11 of http://arxiv.org/abs/hep-ph/0611266
        complex<double> H1 = 9.0*pow(consts::z2,2.0)/10.0
        + 2.0*consts::z3*log(x)
        + consts::z2*chaplin::HPL(0,0,x)
        + 1.0/4.0*chaplin::HPL(0,0,0,0,x)
        + 7.0/2.0*chaplin::HPL(0,1,0,0,x)
        - 2.0*chaplin::HPL(0,-1,0,0,x)
        + 4.0*chaplin::HPL(0,0,-1,0,x)
        - chaplin::HPL(0,0,1,0,x);
        
        complex<double> F2la = 36.0*x/pow(x-1.0,2.0)
        - 4.0*x*(1.0-14.0*x+x*x)/pow(x-1.0,4.0)*consts::z3
        - 4.0*x*(1.0+x)*log(x)/pow(x-1.0,3.0)
        - 8.0*x*(1.0+9.0*x+x*x)/pow(x-1.0,4.0)*chaplin::HPL(0,0,x)
        +2.0*x*(3.0+25.0*x-7.0*x*x+3.0*pow(x,3.0))/pow(x-1.0,5.0)*chaplin::HPL(0,0,0,x)
        + 4.0*x*(1.0+2.0*x+x*x)/pow(x-1.0,4.0)*
            (consts::z2*log(x) + 4.0*chaplin::HPL(0,-1,0,x) - chaplin::HPL(0,1,0,x))
        + 4.0*x*(5.0-6.0*x+5.0*pow(x,2.0))/pow(x-1.0,4.0)*chaplin::HPL(1,0,0,x)
        - 8.0*x*(1.0+x+pow(x,2.0)+pow(x,3.0))/pow(x-1.0,5.0)*H1;
        
        
       
        
        //: putting the above together: we keep here N arbitrary for readability reasons only
        double N=3.0;
        double Cf=(N*N-1.0)/2.0/N;
        complex<double> res =1.0/2.0
                            * (Cf * F2la  +  N * G2lCA_onehalf
                               + Cf * scheme_dependent_coeff * F2lb(x)
                               )
                            / (-2.0/3.0);
        
        return res;
        
    }

    complex<double> born_exact_summed_over_quarks(const CModel& model)
    {
        
        complex<double> ME ;
        
        for (int i=0;i<model.quarks.size();i++)
            {
            ME = ME + model.quarks[i]->Y() * born(model.quarks[i]->X());
            }
        return(ME);
    }

    complex<double> born(complex<double> x)
    {
        
        //: the expression below goes to 1 as mq->infty, i.e. as x->1
        complex<double > res=(-3.0)*x/pow(1.0-x,2.0)*
                    (2.0-pow(1.0+x,2.0)/pow(1.0-x,2.0)*chaplin::HPL(0,0,x));
        
        if (res!=res){
            cout << "nan in Born : x = " << x << endl;
        }
        
        return res;
    }

    double sum_of_abs_sq_of_Aqi(const double &z,const double & lambda,const CModel& model)
    {
        complex<double> ME ;
        double mh = model.higgs.m();
        complex<double> A1(0.0,0.0);
        complex<double> A2(0.0,0.0);
        complex<double> A3(0.0,0.0);
        complex<double> A4(0.0,0.0);
        for (int i=0;i<model.quarks.size();i++)
            {
            double Yq = model.quarks[i]->Y();
            if (Yq!=0.0)
                {
                complex<double> Wq = model.quarks[i]->Wq();
                complex<double> mq_cplx = sqrt(
                                model.quarks[i]->cm_sq());
                complex<double> factor = Yq * Wq * 3.0/32.0 /pow(mh,2.0);
                A1 += Aq1(z,lambda,mq_cplx,mh) * factor;
                A2 += Aq2a(z,lambda,mq_cplx,mh) * factor;
                A3 += Aq2b(z,lambda,mq_cplx,mh) * factor;
                A4 += Aq2c(z,lambda,mq_cplx,mh) * factor;
                //cout<<"\n%%% "<<model.quarks[i]->name<<" "<<z<<" "<<lambda;
                }
            }
            return(    pow(abs(A1),2.0)
                    + pow(abs(A2),2.0)
                    + pow(abs(A3),2.0)
                    + pow(abs(A4),2.0));
    }

    complex<double> sum_of_Aqqgh(const double& y,const CModel& model)
    {
        complex<double> res(0.0,0.0);
        for (int i=0;i<model.quarks.size();i++)
        {
            double Yq = model.quarks[i]->Y();
            if (Yq!=0.0)
            {
                complex<double> Wq = model.quarks[i]->Wq();
                complex<double> factor = Yq * Wq * 3.0 / 4.0;
                res += Aqqgh_cpp(y,model.quarks[i]->X()) * factor;
            }
        }
        return res;
    }
    
    
    double sum_of_abs_sq_of_Aqqgh(const double& y,const CModel& model)
    {
        return pow(abs(sum_of_Aqqgh(y,model)),2.0);
    }
    
    complex<double> Aqqgh_cpp(const double& y,const complex<double>& x)
    {
        // y = -1/zbar/lambdabar in the case of qg->qH
        complex<double> x12 = sqrt(1.0+4.0*x*y/pow(1.0-x,2.0));
        complex<double> TS = spec_triaf(x12);
        complex<double> xmh = sqrt(1.0-(1.0-x12*x12)/y);
        complex<double> TH = spec_triaf(xmh);
        double ybar = 1.0-y;
        complex<double> res = y/ybar * (1.0/3.0 * pow(1.0-x,2.0)/x * born(x)
                                        - 1.0 / ybar * 2.0 * (bub_aux(x12)-bub_aux(xmh))
                                        + 1.0 / ybar * 4.0*x/pow(1.0-x,2.0) * (TS * y -TH)
                                        -  TS);
        return res;
    }
    

    complex<double> Aq1(const double& z, const double& lambda,
                        const complex<double> & m, const double& mh)
    {
        const double s = pow(mh,2.) / z;
        const double t = -s*(1-z)*lambda;
        const double u = -s*(1-z)*(1-lambda);
        
        const double c11 = (-16*(pow(s,3)*t + pow(s,2)*pow(t,2) + s*pow(t,3)
                                 + pow(s,3)*u + pow(s,2)*t*u + s*pow(t,2)*u
                                 + pow(t,3)*u + pow(s,2)*pow(u,2) + s*t*pow(u,2)
                                 + pow(t,2)*pow(u,2) + s*pow(u,3) + t*pow(u,3)
                                 )
                            )/((s + t)*(s + u)*(t + u));
        
        const complex<double> c7 =  (-8*s*t)/u - (8*s*u)/t - (8*t*u)/s
                            - CCtri(m,s,t,u)
                            - CCtri(m,t,u,s)
                            - CCtri(m,u,s,t);
        
        return           c11
                            + boxf(m,s,t,u)*CCbox(s,t,u)
                            + boxf(m,t,u,s)*CCbox(t,u,s)
                            + boxf(m,u,s,t)*CCbox(u,s,t)
                            + bubf(m,s,s + t + u)*CCbub(s,t,u)
                            + bubf(m,t,s + t + u)*CCbub(t,u,s)
                            + bubf(m,u,s + t + u)*CCbub(u,s,t)
                            + CCtri(m,s,t,u)*triaf(m,s)
                            + CCtri(m,t,u,s)*triaf(m,t)
                            + CCtri(m,u,s,t)*triaf(m,u)
                            + c7*triaf(m,s + t + u);

    }


    complex<double> Aq2a(const double& z,const double& lambda,
                         const complex<double>& m,const double& mh)
    {
        const double s = pow(mh,2.) / z;
        const double t = -s*(1-z)*lambda;
        const double u = -s*(1-z)*(1-lambda);
        return Aq2(s,t,u,m);
    }

    complex<double> Aq2b(const double& z,const double& lambda,
                         const complex<double>& m,const double& mh)
    {
        const double s = pow(mh,2.) / z;
        const double t = -s*(1-z)*lambda;
        const double u = -s*(1-z)*(1-lambda);
        return Aq2(t,u,s,m);
    }

    complex<double> Aq2c(const double& z,const double& lambda,
                         const complex<double>& m,const double& mh)
    {
        const double s = pow(mh,2.) / z;
        const double t = -s*(1-z)*lambda;
        const double u = -s*(1-z)*(1-lambda);
        return Aq2(u,s,t,m);
    }


    complex<double> Aq2(const double& s,const double& t,const double& u,const complex<double>& m)
    {
        const double u_sq = u*u;
        const double t_sq = t*t;
        const double s_sq = s*s;
        const double u_cube = u*u*u;
        const double s_cube = s*s*s;
        
        const complex<double> cg1 = (-16. * pow(m,2.) - 2. * (-u_sq - t * u + 2. * s * t) / u);
        const complex<double> cg2 = (-16. * pow(m,2.) + 2. * (s * t + s_sq - 2. * t * u) / s);
        const complex<double> cg3 =  (16. * pow(m,2.) + 4. * s * u / t);
        const complex<double> cg4 = 8. * s * (-u + t) / (u + t);
        const complex<double> cg5 = 8. * (u_sq + 4. * s * u + s_sq) / pow(u + s,2.) * t;
        const complex<double> cg6 = 8. * (s - t) / (t + s) * u;
        const complex<double> cg7 =  (   -16. * (8. * s * t_sq * u + 3. * t_sq * u_sq
              + 5. * s_sq * u * t + u_cube * t + t * s_cube
              - s * u_cube + 3. * s_sq * t_sq - s_cube * u
              + 5. * s * u_sq * t
              ) / pow(u + s,2.) / (t + s) / (u + t)* pow(m,2.)
          - 4. * (-s_sq * t_sq * u + 2. * s_cube * t_sq
                  - s * t_sq * u_sq + 2. * t_sq * u_cube
                  - 2. * t * u_sq * s_sq - 2. *s_sq * u_cube
                  - 2. * s_cube * u_sq
                   ) / (u + s) / s / t / u);
        const complex<double> cg8 = (16. * (-u + t) / (u + t) * pow(m,2.)
                  + 4. * (-t_sq * u + 2. * s * t_sq
                      - 2. * u_sq * s)
                               / t / u);
        const complex<double> cg9 = (16. * (u_sq + 4. * s * u + s_sq) / pow(u + s,2.)
                  * pow(m,2.)
              + 4. * (2. * t * s_cube - 2. * s_sq * u_sq
                      + s_sq * u * t + s * u_sq * t
                      + 2. * u_cube * t)
                                      / u / s / (u + s));
                                      
        const complex<double> cg10 = (-16. * (s - t) / (t + s) * pow(m,2.)
              - 4. * (2. * u * s_sq + s * t_sq
                       - 2. *t_sq * u) / s / t);
        const complex<double> cg11 = 16. * (s * u + s * t + t * u) / (u + s);
                                              
                                              
        return  cg1 * boxf(m, s, t, u)
            +  cg2 * boxf(m, t, u, s)
            +  cg3 * boxf(m, u, s,t)
            +  cg4 * bubf(m, s, s + t + u)
            +  cg5 * bubf(m, t, s +t + u)
            -  cg6 * bubf(m, u, s + t + u)
            +  cg7 * triaf(m, s + t + u)
            +  cg8 * triaf(m, s)
            +  cg9 * triaf(m, t)
            + cg10 * triaf(m, u)
            - cg11;
    }


    // coefficients of Aq1
    double CCbox(const double& s, const double& t, const double& u)
    {
        return (2*(2*s*t + s*u + t*u))/u;
    }

    complex<double> CCtri(const complex<double>& m, const double& s,
                          const double& t, const double& u)
    {
        return (16.*pow(m,2)*(pow(t,2) + pow(u,2)))/pow(t + u,2)
        - (4.*(2.*s*pow(t,3)+ 3.*s*pow(t,2)*u + pow(t,3)*u + 3.*s*t*pow(u,2)
               + 2.*s*pow(u,3) + t*pow(u,3))
           )/(t*u*(t + u));
    }

    double CCbub(const double& s, const double& t, const double& u)
    {
        return (-8*s*(pow(t,2) + 4*t*u + pow(u,2)))/pow(t + u,2);
    }


    // master integrals
    complex<double> bubf(const complex<double>& m, const double& s1, const double& s2)
    {
        return bubble(m,s1)-bubble(m,s2);
    }

    complex<double> bubble(const complex<double>& m, const double& s)
    {
        complex<double> x = sqrt(1.-4.*m*m/s);
        return bub_aux(x);
    }

    complex<double> bub_aux(const complex<double>& x)
    {
        complex<double> Ie(0.0,1e-15);
        return x*log((x-1.)/(x+1.)+Ie);
    }

    // triangle: triaf(m,s) = C_0(m,s)*s
    complex<double> triaf(const complex<double>& m, const double& s)
    {
        complex<double> x = sqrt(1.-4.*m*m/s);
        complex<double> Ie(0.0,1e-15);
        return 0.5 * pow(log((x-1.)/(x+1.)+Ie),2.);
    }
    // same triangle but written as a function of one variable, 'x'
    complex<double> spec_triaf(const complex<double>& x)
    {
        complex<double> Ie(0.0,1e-15);
        return 0.5*pow(log( (x-1.0)/(x+1.0) + Ie ),2.0);
    }




    



    //------ Box according to Glover & Baur, Nucl.Phys.B339(1990)
    //       Note that boxf(m,s,t,u) = D0(m,s,t,u)*s*t
    //       Checks whether the mass is real (no complex mass scheme) -> Glover & Baur formula
    //       or we are in the complex mass scheme -> avd_olo_d0c
    complex<double> boxf(const complex<double>& m,
                         const double& s,
                         const double& t,
                         const double& u)
    {
        complex<double> res;
        if (abs(imag(m))<1e-10)
        {
        double Q = s+t+u;

         res =    - Box_Intv2(s,t,u,m,Q)
             + Box_Intv2(s,t,u,m,s)
             + Box_Intv2(s,t,u,m,t);
        }
        else
        {
            double Q = s+t+u;
            
            complex<double> andre(0.0,0.0);
            complex<double> zero(0.0,0.0);
            complex<double> cQ(Q,0.);
            complex<double> cs(s,0.);
            complex<double> ct(t,0.);
            
            
            res = h_exact::boxf_complex_masses(m,s,t,u) * s * t;

        }
        return res;
        
    }


    // -----  The basic integral of eq.B9 in Glover & Baur, Nucl.Phys.B339(1990)
    //       version 2 : all imaginary parts are set explicitly by an Ie
    //       All results are cross checked for s*t*u>0
    complex<double> Box_Intv2(const double& s,
                              const double& t,
                              const double & u,
                              const complex<double> & m,
                              const double& Q)
    {
        if (s*t*u<0.0)
        {
            cout<<"[Box_Intv2]: u/s/t < 0 - risky kinematical region "
                <<"for the box implementation. Check against other libraries."
                <<endl;
        }
        const complex<double> Ie(0.0,1e-15);
        const complex<double> r = 1. + 4.*(m*m-Ie)*u/t/s;
        const complex<double> x_tau = sqrt(r);
        

        const complex<double> ry=1.-4.*(pow(m,2.)-Ie)/Q;
        const complex<double> y_tau = sqrt(ry);
        const complex<double> x_plus = 0.5*(1.+x_tau);
        const complex<double> x_minus = 0.5*(1.-x_tau);
        const complex<double> y = 0.5*(1. + y_tau);
                             
        const complex<double> the_log_1 = log((x_tau-1.)/(x_tau+1.));
                             
        const complex<double> the_log_2 = log(1. + (Q-Ie)*u/s/t);
        
        
        return 2./x_tau * (
                         chaplin::HPL(0,1,x_minus/(x_minus - y )-Ie)
                        -chaplin::HPL(0,1,x_plus /(x_plus  - y )+Ie)
                        +chaplin::HPL(0,1,x_minus/( y - x_plus )+Ie)
                        -chaplin::HPL(0,1,x_plus /( y -x_minus )-Ie)
                         +the_log_1*the_log_2
                                );
                               

    }



//: new implementation of the special case box we need for nlo QCD corrections to gluon fusion with exact quark masses that are COMPLEX.

    double theta(const double& x) {
        if ( x >= 0.0 ) return 1.0;
        else return 0.0;
    }

    complex<double> nn(const complex<double>& a, const complex<double>& b) {
        double Ima = a.imag();
        double Imb = b.imag();
        const complex<double> ab = a * b;
        double Imab = ab.imag();
        double one_or_zero = theta(-Ima)*theta(-Imb) * theta(Imab)
                            -theta(Ima) * theta(Imb) * theta(-Imab);
        return complex<double> (0.0,2.*consts::Pi * one_or_zero );
    }

    complex<double> RR(const complex<double> & y, const complex<double> & z) {
        const complex<double> i_epsilon (0.0,1e-15);
        return log(y-z) * log( (1.-y) / (-y) -i_epsilon )
            + chaplin::HPL(0,1, (y)/(y-z) )
            - chaplin::HPL(0,1, (y-1.)/(y-z) )
            + nn(-z,1./(y-z)) * log( y/(y-z) )
            - nn(1.-z,1./(y-z)) * log( (y-1.)/(y-z) ) ;
    }

    complex<double> SS(const double& a, const complex<double>& m_sq, const complex<double>& y) {
        const complex<double> z1 = (1. + sqrt(1. - 4. * m_sq / a) ) / 2.;
        const complex<double> z2 = (1. - sqrt(1. - 4. * m_sq / a) ) / 2.;
        const complex<double> i_epsilon (0.0,1e-15);
        return ( log( a - i_epsilon ) + nn( -z1, -z2) ) * log( (y - 1.)/y )
                + RR(y,z1) + RR(y,z2);
        
    }

    complex<double> boxf_complex_masses(const complex<double>& m,
                                        const double& s,
                                        const double& t,
                                        const double& u) {
        if (abs(m.imag()) < 1e-16) {
            cout << " Error: you called the wrong implementation for the box: "
                << " boxf_complex_masses assumes the mass has non-zero imaginary part"
                    <<endl;
            exit(EXIT_FAILURE);
        }
        double Qsq = s + t + u ;
        complex<double> m_sq = m*m;
        complex<double> Delta = 1. + 4. * m_sq * u / s / t;
        complex<double> y1 = (1. - sqrt(Delta))/2.;
        complex<double> y2 = (1. + sqrt(Delta))/2.;
        const complex<double> i_epsilon (0.0,1e-15);
        return 1. / s / t / sqrt(Delta) * (
                            SS(s,m_sq,y1) - SS(s,m_sq,y2)
                            + 2. * log(m_sq)
                                 * log( (sqrt(Delta) -1.) / (sqrt(Delta) + 1.) )
                            + SS(t,m_sq,y1) - SS(t,m_sq,y2)
                            - SS(Qsq,m_sq,y1) + SS(Qsq,m_sq,y2)
                            );
        
    }

}






