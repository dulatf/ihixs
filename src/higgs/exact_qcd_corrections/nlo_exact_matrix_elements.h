#ifndef EHIXS_GLUON_FUSION_NLOEXACTMATRIXELEMENTS_H
#define EHIXS_GLUON_FUSION_NLOEXACTMATRIXELEMENTS_H

#include <complex>
using namespace std;
class CModel;

namespace h_exact {


    class ExactRegularCalculator{
    public:
        ExactRegularCalculator(CModel* model){_model=model;SetUp();}
        void SetUp();
        double gg_reg_Log0(const double& z, const double& lambda);
        double gg_reg_Log0_hard(const double& z, const double& lambda);
        double gg_reg_Log0_collinear(const double& z, const double& lambda);
        double gg_reg_Log1(const double& z, const double& lambda);
        double qg_reg_Log0(const double& z,const double& lambda);
        double qg_reg_Log0_hard(const double& z,const double& lambda);
        double qg_reg_Log0_collinear(const double& z,const double& lambda);
        double qg_reg_Log1(const double& z,const double& lambda);
        double qqb_reg_Log0(const double& z,const double& lambda);
        
        void compute_soft();
        double born_squared(){return _B_sq;}
        double soft(){return _soft;}
    private:
        CModel* _model;
        double _mh;
        vector< complex<double> > _factor;
        vector< complex<double> > _factorqqgh;
        vector< complex<double> > _mq_cplx;
        vector< complex<double> > _Xfactor;
        double _B_sq;
        double _soft;
    private:
        double sum_of_abs_sq_of_Aqi(const double &z,const double & lambda);
        complex<double> sum_of_Aqqgh(const double& y);
        double sum_of_abs_sq_of_Aqqgh(const double& y);

        
    };
    
    complex<double> ggf_exact_virtual_ep0(const complex<double> & x,
        const complex<double>& scheme_dependent_coeff);
    double gg_reg_exact_nlo(const double& z, const double& lambda,
                            const double& L,const CModel& model);
    
    double gg_reg_exact_nlo_L0(const double& z, const double& lambda,
                            const CModel& model);
    double gg_reg_exact_nlo_L0_hard_only(const double& z,
                                         const double& lambda,
                                         const CModel& model);
    double gg_reg_exact_nlo_L1(const double& z, const double& lambda,
                               const CModel& model);
    
    double qg_reg_exact_nlo(const double& z, const double& lambda,
                            const double& L,const CModel& model);
    double qg_reg_exact_nlo_L0(const double& z, const double& lambda,
                            const CModel& model);
    double qg_reg_exact_nlo_L0_hard_only(const double& z, const double& lambda,
                               const CModel& model);
    double qg_reg_exact_nlo_L1(const double& z, const double& lambda,
                            const CModel& model);
    
    double qqb_reg_exact_nlo(const double& z, const double& lambda,
                            const double& L,const CModel& model);
    double qqb_reg_exact_nlo_L0(const double& z, const double& lambda,
                             const CModel& model);
   
    
    complex<double> F2lb(const complex<double> & x);




    typedef complex<double> (*ptr_to_Aq)(const double& z,
                                    const double& lambda,
                                    const complex<double>& M,
                                    const double& QQQ);

    double abs_sq_of_sum_over_quarks_of(ptr_to_Aq,const double & z,
                                    const double & lambda,const CModel& model);
    double sum_of_abs_sq_of_Aqi(const double &z,const double & lambda,
                                const CModel& model);
    double sum_of_abs_sq_of_Aqqgh(const double& z,const CModel& model);
    complex<double> sum_of_Aqqgh(const double& y,const CModel& model);
    
    complex<double> Aqqgh_cpp(const double& z,
                              const complex<double>&x );
    
    
    complex<double> born_exact_summed_over_quarks(const CModel& model);
    complex<double> born(complex<double> x);


    complex<double> Aq1(const double& z, const double& lambda,
                        const complex<double> & m, const double& mh);
    complex<double> Aq2a(const double& z,const double& lambda,
                         const complex<double>& m,const double& mh);
    complex<double> Aq2b(const double& z,const double& lambda,
                         const complex<double>& m,const double& mh);
    complex<double> Aq2c(const double& z,const double& lambda,
                         const complex<double>& m,const double& mh);
    complex<double> Aq2(const double& s,const double& t,const double& u,const complex<double>& m);

    // coefficients of Aq1
    double CCbox(const double& s, const double& t, const double& u);
    complex<double> CCtri(const complex<double>& m, const double& s,
                          const double& t, const double& u);
    double CCbub(const double& s, const double& t, const double& u);

    // masters
    complex<double> boxf(const complex<double>& m,const double& s,const double& t,const double& u);
    complex<double> triaf(const complex<double>& m, const double& s);
    complex<double> spec_triaf(const complex<double>& x);

    complex<double> bubf(const complex<double>& m, const double& s1, const double& s2);
    complex<double> bubble(const complex<double>& m, const double& s);
    complex<double> bub_aux(const complex<double>& x);
    
    complex<double> Box_Intv2(const double& s,const double& t,const double & u,
                              const complex<double> & m,const double& Q);
    //: new implementation
    complex<double> boxf_complex_masses(const complex<double>& m,
                                        const double& s,
                                        const double& t,
                                        const double& u);
    complex<double> SS(const double& a, const complex<double>& m_sq, const complex<double>& y);

    
    
    
}
#endif
