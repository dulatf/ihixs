/** testing UserInterface.*:
 *
 * Achilleas Lazopoulos, lazopoli@phys.ethz.ch
 */

#include <iostream>
#include <cmath>



using namespace std;

#include "gtest/gtest.h"
#include "psi.h"

PsiFinder *PsiFinder::s_instance = 0;



double diff(const complex<double> res, const complex<double> exp){
    
    if (res==exp) return 0.0;
    
    complex<double> diff = res-exp;
    
    if (abs(exp)<1e-15) return abs(diff);
    else return abs(diff)/abs(exp);
    
}

TEST(PsiTest,a1)
{
    complex<double> z(2.5,0.343426);
    complex<double> res = PsiFinder::psi()->get(z);
    complex<double> exp = complex<double>(0.716907 , 0.166913 );
    double err=1e-5;
    
    
    
    EXPECT_LT(abs(real(res)-real(exp))/abs(real(exp)),err)<<" Mismatch in real part: res="<<res<<"| expecting="<<exp<<"|"<<endl;
    EXPECT_LT(abs(imag(res)-imag(exp))/abs(imag(exp)),err)<<" Mismatch in imag part: res="<<res<<"| expecting="<<exp<<"|"<<endl;
}

TEST(PsiTest,a2)
{
    complex<double> z(3.5,0.343426);
    complex<double> res = PsiFinder::psi()->get(z);
    complex<double> exp = complex<double>(1.1095 , 0.112982 );
    double err=1e-5;
    
    
    
    EXPECT_LT(abs(real(res)-real(exp))/abs(real(exp)),err)<<" Mismatch in real part: res="<<res<<"| expecting="<<exp<<"|"<<endl;
    EXPECT_LT(abs(imag(res)-imag(exp))/abs(imag(exp)),err)<<" Mismatch in imag part: res="<<res<<"| expecting="<<exp<<"|"<<endl;
}

TEST(PsiTest,a3)
{
    complex<double> z(4.5,0.343426);
    complex<double> res = PsiFinder::psi()->get(z);
    complex<double> exp = complex<double>(1.39249 , 0.085215 );
    double err=1e-5;
    
    
    
    EXPECT_LT(abs(real(res)-real(exp))/abs(real(exp)),err)<<" Mismatch in real part: res="<<res<<"| expecting="<<exp<<"|"<<endl;
    EXPECT_LT(abs(imag(res)-imag(exp))/abs(imag(exp)),err)<<" Mismatch in imag part: res="<<res<<"| expecting="<<exp<<"|"<<endl;
}


TEST(PsiTest,b1)
{
    complex<double> z(2.5,124.983);
    complex<double> res = PsiFinder::psi()->get(z);
    complex<double> exp = complex<double>(4.8283 , 1.5548  );
    double err=1e-5;
    
    
    
    EXPECT_LT(abs(real(res)-real(exp))/abs(real(exp)),err)<<" Mismatch in real part: res="<<res<<"| expecting="<<exp<<"|"<<endl;
    EXPECT_LT(abs(imag(res)-imag(exp))/abs(imag(exp)),err)<<" Mismatch in imag part: res="<<res<<"| expecting="<<exp<<"|"<<endl;
}

TEST(PsiTest,b2)
{
    complex<double> z(3.5,124.983);
    complex<double> res = PsiFinder::psi()->get(z);
    complex<double> exp = complex<double>(4.82846 , 1.5468  );
    double err=1e-5;
    
    
    
    EXPECT_LT(abs(real(res)-real(exp))/abs(real(exp)),err)<<" Mismatch in real part: res="<<res<<"| expecting="<<exp<<"|"<<endl;
    EXPECT_LT(abs(imag(res)-imag(exp))/abs(imag(exp)),err)<<" Mismatch in imag part: res="<<res<<"| expecting="<<exp<<"|"<<endl;
}

TEST(PsiTest,b3)
{
    complex<double> z(4.5,124.983);
    complex<double> res = PsiFinder::psi()->get(z);
    complex<double> exp = complex<double>(4.82869 ,1.5388 );
    double err=1e-5;
    
    
    
    EXPECT_LT(abs(real(res)-real(exp))/abs(real(exp)),err)<<" Mismatch in real part: res="<<res<<"| expecting="<<exp<<"|"<<endl;
    EXPECT_LT(abs(imag(res)-imag(exp))/abs(imag(exp)),err)<<" Mismatch in imag part: res="<<res<<"| expecting="<<exp<<"|"<<endl;
}



TEST(PsiTest,Newb3)
{
    complex<double> z(17.5,124.983);
    complex<double> exp = PsiFinder::psi()->get(z);
    complex<double> res = psi(z);
    double err=1e-5;
    
    
    
    EXPECT_LT(abs(real(res)-real(exp))/abs(real(exp)),err)<<" Mismatch in real part: res="<<res<<"| expecting="<<exp<<"|"<<endl;
    EXPECT_LT(abs(imag(res)-imag(exp))/abs(imag(exp)),err)<<" Mismatch in imag part: res="<<res<<"| expecting="<<exp<<"|"<<endl;
}


TEST(PsiTest,Kolbig1)
{
    complex<double> z(-13.0,2.0);
    complex<double> exp(2.613758858614923,2.994600955642857);
    complex<double> res = psi(z);
    double err=1e-15;
    
    EXPECT_LT(abs(res-exp),err) << " res = " << res
    << " exp = " << exp << endl;
}

TEST(PsiTest,Kolbig2)
{
    complex<double> z(0.0,8.0);
    complex<double> exp(2.080745674911801,1.633296326794897);
    complex<double> res = psi(z);
    double err=1e-15;
    
    EXPECT_LT(abs(res-exp),err) << " res = " << res
    << " exp = " << exp << endl;
}

TEST(PsiTest,Kolbig3)
{
    complex<double> z(3.0,0.0);
    complex<double> exp(0.922784335098467,0.0);
    complex<double> res = psi(z);
    double err=1e-15;
    
    EXPECT_LT(abs(res-exp),err) << " res = " << res
    << " exp = " << exp << endl;
}


TEST(PsiTest,Kolbig4)
{
    complex<double> z(15.0,-5.0);
    complex<double> exp( 2.7304638296862949658807696179355,
                        -0.3319504266337825064013487331467);
    complex<double> res = psi(z);
    //complex<double> res = PsiFinder::psi()->get(z);
    double err=1e-15;
    
    EXPECT_LT(abs(res-exp),err) << " res = " << res
                                << " exp = " << exp << endl;

}


class PsiMultiTest: public ::testing::Test {
protected:
    virtual void SetUp()
    {
        zero = complex<double> (0.0,0.0);
        one  = complex<double> (1.0,0.0);
        minus_one =  complex<double> (-1.0,0.0);
        zvals = vector<complex<double> >{ complex<double>(0.5,0.0),
            complex<double>(2.,0.0),
            //complex<double>(-1.,0.0),
            //complex<double>(0.0,0.0),
            complex<double>(0.5,0.5),
            complex<double>(2.0,2.0),
            complex<double>(0.02,0.003),
            complex<double>(0.2,0.1),
            complex<double>(-0.2,0.0),
            complex<double>(-10.2,0.0),
            complex<double>(10.34,0.0),
            complex<double>(1.0,1e-6),
            complex<double>(1.0,-1e-6),
            complex<double>(1.0,1e-3),
            complex<double>(1.0,-1e-3),
            complex<double>(0.8,0.1),
            complex<double>(0.01,0.1)
        };
        int N=10;
        for (int i = 0; i < N; i++){
            for (int j = 0; j < N; j++){
                zvals.push_back(complex<double>(-10.123 + double(i)/N*20.01,
                                                -10.456 + double(j)/N*20.03)
                                );
            }
        }
        
        err=1e-14;
        
    }
    
    vector<complex<double> >zvals;
    complex<double> zero,one,minus_one;
    double err;
};


TEST_F(PsiMultiTest,all)
{
   for (auto z : zvals) {
        complex<double> res = psi(z);
        complex<double> exp = PsiFinder::psi()->get(z);
        EXPECT_LT(diff(res , exp ),err)
        << " Mismatch for z= " << z
        <<" : res=" << res
        << "| expecting=" << exp << "|" << endl;
        
    }
    
    
}
int main(int argc, char**argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return  RUN_ALL_TESTS();
    return 0;
}




















