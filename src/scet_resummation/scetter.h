#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <vector>
#include <string>
#include <iterator>
#include <sys/types.h>
#include <sys/stat.h>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <complex>

#include <boost/math/special_functions/polygamma.hpp>
#include <boost/math/special_functions/gamma.hpp>

#include "LHAPDF/GridPDF.h"
#include "cuba.h"

using namespace std;

enum PDFs{
    MSTW08=0,
    ABM=1,
    NNPDF=2,
    PDF4LHC=3
};


class SCETTER {
    
public:
    
    vector<double> beta,gammaC,gammaS,gammag,gammaABNY;
    vector<double> GlobalH,GlobalU,GlobalW;
    double Pi,EulerGamma;
    double Gf;
    double IsVegas;
    
    double globaltau,globalmh,globalmt,globalmut,globalmuh,globalx,globalorder;
    double globalmode,globaldist,globalmus,globalmuf,globalpisq,globalLt;
    
    double gammaC3fudge;
    double gammaG3fudge;
    double gammaS3fudge;
    int alpha_order;
    double Lambda_QCD;
    
    double globalL1;
    double globalL1p;
    double globalL1pp;
    
    
    LHAPDF::PDF* pdf;
    double alphaS;
    
    
    SCETTER(){
        pdf=NULL;
        gammaC3fudge=0.0;
        gammaG3fudge=0.0;
        gammaS3fudge=0.0;
        alpha_order=3;
        Lambda_QCD=2.50008;
        
        IsVegas=0;
        Gf=0;
        
        
        globalL1=0;
        globalL1p=0;
        globalL1pp=0;
        
        
        Pi=3.14159265358979323846264338328;
        EulerGamma=0.577215664901532860606512090082;
        
        //BetaFunction
        beta.push_back(1.91666666666666666666666666667);
        beta.push_back(2.41666666666666666666666666667);
        beta.push_back(2.82667824074074074074074074074);
        beta.push_back(18.8521731593394402719352219415);
        
        
        //Cusp anomalous dimension
        gammaC.push_back(1.0);
        gammaC.push_back(1.72704334417210478973582169448);
        gammaC.push_back(2.80321913905030613872095606875);
        gammaC.push_back(gammaC3fudge);
        
        //Gluon anomalous dimension
        gammag.push_back(1.91666666666666666666666666667);
        gammag.push_back(10.7805507629939280931148992569);
        gammag.push_back(30.3660101909229956040514260678);
        gammag.push_back(gammaG3fudge);
        
        
        //Soft anomalous dimension
        gammaS.push_back(0.0);
        gammaS.push_back(-1.17729763681337296410350664360);
        gammaS.push_back(-20.6033259705267120452457226618);
        gammaS.push_back(gammaS3fudge);
        
        
        ///WTF ABNY
        gammaABNY.push_back(gammaS[0]-2.0*beta[0]);
        gammaABNY.push_back(gammaS[1]-4.0*beta[1]);
        gammaABNY.push_back(gammaS[2]-6.0*beta[2]);
        gammaABNY.push_back(gammaS[3]-8.0*beta[3]);
        
        
    };
    
    //PDF Stuff
    double GetAlpha(double Q);
    double GetValue(double x, double Q,int parton);
    double GetDerivative(double x, double Q,int parton);
    double GetDDerivative(double x, double Q,int parton);
    double GetDDDerivative(double x, double Q,int parton);
    void InitiatePDF(const string& pdfset);
    
    
    //Random Functions
    double Zeta(int i);
    vector<double> Multiply(vector<double> a,vector<double>b);

    
    //AlphaS Stuff
    complex<double> alpha_res_complex(double Lambda_QCD,double scale,int order);
    double alpha_res_real(double Lambda_QCD,double scale,int order);
    void SetLambdaQCD(double a0, double mu0,int order);
    void SetAlphaOrder(int order, double alpha0, double mu0);
    complex<double> alpha(double mu1,int order,bool pisq);
    
    //Anomalous dimensions
    vector<complex<double> > Ag(double mu0,double mu1,vector<double> gamma,bool pisq);
    vector<complex<double> > Sg(double mu0,double mu1,vector<double> gamma,bool pisq);
    
    //Hard Function
    vector<double> Hard(double mH,double mu1,bool pisq);
    
    //Transfer Function
    vector<double> U(double mh,double muf,double mus,double muh,double mut,double mt,bool pisq);
    
    //Wilson Coefficient
    vector<double> WilsonCoeff(double mt,double mut);
    
    //Soft Function Stuff
    //Fit function for the soft scale
    double musfit(double x)
    {
        return -45.51809781055389 - 0.30889436395050324*x + 0.0004072448965095652*pow(x,2) + 28.54248837013623*log(x) + 0.7213561312695667*pow(log(x),2);
    };
    vector<double> GetLD(vector<double> DD,double L1,double L1p,double L1pp,double L1ppp,double x,double tau,double Lt,int mode);
    vector<double> Soft(double mu1,vector<double> LD);
    double GetX(double mus,double muf,int order);
    
    //Stuff that is needed to integrate numerically + distributions
    double Distribution(double z,double mh,double mus,double x,double Lz,double L1,double L1p,double L1pp,double L1ppp,int mode,int num);
    double func(double x,double lowbound);
    double Lumi(const double xx[]);
    double LumiP(const double xx[]);
    double LumiPP(const double xx[]);
    double Integrand(const double xx[]);
    vector<double> IntegrateXS();
    
    
    //XS
    double GetBorn(double mh);
    vector<double> CreateXS(double mh, double mt,double muf,double muh,double mus,double mut,double tau,bool pisq, double sqrts);
    vector<double> DeltaSCET(double mh, double mt,double GF,double muf,double tau);

};

