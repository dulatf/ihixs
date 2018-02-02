/** testing UserInterface.*:
 *
 * Achilleas Lazopoulos, lazopoli@phys.ethz.ch
 */

#include <iostream>
#include <cmath>



using namespace std;

#include "gtest/gtest.h"
#include "higgs_eft.h"
//#include "inclusive_process.h"
#include "math.h"
class HiggsEFTMe: public ::testing::Test {
protected:
    virtual void SetUp()
        {
            err=1e-4;
            z=0.72021484375;
            lz = log(1.-z);
        }
    double err;
    double z;
    double lz;
};


// testing the log(muf/mh) *INDEPENDENT* terms of delta and plus
// partonic cross sections within the Higgs eft, aka Egg_LO, NLO, NNLO and N3LO
// against ihixs3
//
// the delta pieces are pure numbers
// the plus pieces are pure numbers times log(1-z)^a
// we have tested here with z=0.72021484375

TEST_F(HiggsEFTMe,LO_delta)
{
    const double res = HEFT::n_LO_delta();
    const double exp = 1.0;
    EXPECT_LT(abs((res-exp)/exp),err);
}

TEST_F(HiggsEFTMe,NLO_delta)
{
    const double res = HEFT::n_NLO_delta();
    const double exp = 9.8696;
    EXPECT_LT(abs((res-exp)/exp),err);
}

TEST_F(HiggsEFTMe,NNLO_delta)
{
    const double res = HEFT::n_NNLO_delta();
    const double exp = 13.61056;
    EXPECT_LT(abs((res-exp)/exp),err);
}

TEST_F(HiggsEFTMe,N3LO_delta)
{
    const double res = HEFT::n_N3LO_delta();
    const double exp = 1124.308;
    EXPECT_LT(abs((res-exp)/exp),err);
}

TEST_F(HiggsEFTMe,N3LO_delta_Logs_muf)
{
    const double L=log(0.34)*2.;
    const double res = HEFT::n_N3LO_delta_L()*L
                    +HEFT::n_N3LO_delta_L2()*L*L
                    +HEFT::n_N3LO_delta_L3()*L*L*L;
    const double exp = -1354.51;//-1379.1391893802765 ;
    cout << "correct = " << res << endl;
    EXPECT_LT(abs((res-exp)/exp),err);
}

TEST_F(HiggsEFTMe,NLO_plus)
{
    const double res = HEFT::n_NLO_D1()*lz;
    const double exp = -15.2848;
    EXPECT_LT(abs((res-exp)/exp),err);
}

TEST_F(HiggsEFTMe,NNLO_plus)
{
    const double res =   HEFT::n_NNLO_D0()
                        +HEFT::n_NNLO_D1()*lz
                        +HEFT::n_NNLO_D2()*lz*lz
                        +HEFT::n_NNLO_D3()*lz*lz*lz;
    const double exp = 161.259;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST_F(HiggsEFTMe,N3LO_plus)
{
    const double res =   HEFT::n_N3LO_D0()
                        +HEFT::n_N3LO_D1()*lz
                        +HEFT::n_N3LO_D2()*lz*lz
                        +HEFT::n_N3LO_D3()*lz*lz*lz
                        +HEFT::n_N3LO_D4()*pow(lz,4.)
                        +HEFT::n_N3LO_D5()*pow(lz,5.);
    
    const double exp = 23173.424870322367;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(QuarkGluonMe,NLO_real)
{
    const double err=1e-4;
    const double z=  0.85171737512117740      ;
    const double L=  -2.4079456086518722      ;
    const double exp = -0.45699212489258401      ;
    const double res = HEFT::qg_nlo_r_lz0(z, L)
                    +HEFT::qg_nlo_r_lz1(z, L);
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;

}

TEST(QuarkGluonMe,NNLO_real)
{
    const double err=1e-4;
    const double z=   4.0067637905890435E-002 ;
    const double L=  -2.4079456086518722      ;
    const double exp =   1956.8745797159809      ;
    const double res = HEFT::qg_nnlo_r_lz0_const(z, L)
                        +HEFT::qg_nnlo_r_lz0_logz(z, L)
                        +HEFT::qg_nnlo_r_lz0_logz_sq(z, L)
                        +HEFT::qg_nnlo_r_lz0_logz_cube(z, L)
                        +HEFT::qg_nnlo_r_lz1(z, L)
                        +HEFT::qg_nnlo_r_lz2(z, L)
                        +HEFT::qg_nnlo_r_lz3(z, L)
                        ;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
    
}
/*
TEST(QuarkGluonMe,DISABLED_N3LO_real_pure_L)
{
    const double err=1e-4;
    const double z=  0.947778320312500000      ;
    const double L=  -2.4079456086518722      ;
    const double exp =  -110.31958404031639      ;
    const double res =
     HEFT::qg_n3lo_r_lz0(z, L)-HEFT::qg_n3lo_r_lz0(z, 0.0)
    +HEFT::qg_n3lo_r_lz1(z, L)-HEFT::qg_n3lo_r_lz1(z, 0.0)
    +HEFT::qg_n3lo_r_lz2(z, L)-HEFT::qg_n3lo_r_lz2(z, 0.0)
    +HEFT::qg_n3lo_r_lz3(z, L)-HEFT::qg_n3lo_r_lz3(z, 0.0)
    +HEFT::qg_n3lo_r_lz4(z, L)-HEFT::qg_n3lo_r_lz4(z, 0.0)
    +HEFT::qg_n3lo_r_lz5(z, L)-HEFT::qg_n3lo_r_lz5(z, 0.0)
    ;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
    
}
*/


TEST(QuarkAntiQuarkMe,NLO_real)
{
    const double err=1e-14;
    const double z=  0.26902993725957480      ;
    const double L=  -2.4079456086518722      ;
    const double exp =   1.7206176577339973      ;
    const double res =HEFT::qqb_nlo_r_lz0(z, L);

    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
    
}

TEST(QuarkAntiQuarkMe,NNLO_real)
{
    const double err=1e-11;
    const double z=  0.38549804687500000      ;
    const double L=  -2.4079456086518722      ;
    const double exp = -0.63860952025809681      ;
    const double res =HEFT::qqb_nnlo_r_lz0_const(z, L)
    +HEFT::qqb_nnlo_r_lz0_logz(z, L)
    +HEFT::qqb_nnlo_r_lz0_logz_sq(z, L)
    +HEFT::qqb_nnlo_r_lz0_logz_cube(z, L)
    +HEFT::qqb_nnlo_r_lz1(z, L)
    +HEFT::qqb_nnlo_r_lz2(z, L)
    ;
    
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
    
}

TEST(QuarkQuarkMe,NNLO_real)
{
    const double err=1e-10;
    const double z=  0.99284246852927549      ;
    const double L=  -2.4079456086518722      ;
    const double exp =  0.22253443807349785      ;
    const double res =HEFT::qq_nnlo_r_lz0_const(z, L)
    +HEFT::qq_nnlo_r_lz0_logz(z, L)
    +HEFT::qq_nnlo_r_lz0_logz_sq(z, L)
    +HEFT::qq_nnlo_r_lz0_logz_cube(z, L)
    +HEFT::qq_nnlo_r_lz1(z, L)
    +HEFT::qq_nnlo_r_lz2(z, L)
    ;
    
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
    
}
/*
TEST(QuarkQuarkPrimeMe,NNLO_real)
{
    const double err=1e-10;
    const double z=   2.1448782891849743E-002 ;
    const double L=  -2.4079456086518722      ;
    const double exp =   2367.8814723140094      ;
    const double res =HEFT::q1q2_nnlo_r_lz0_const(z, L)
    +HEFT::q1q2_nnlo_r_lz0_logz(z, L)
    +HEFT::q1q2_nnlo_r_lz0_logz_sq(z, L)
    +HEFT::q1q2_nnlo_r_lz0_logz_cube(z, L)
    +HEFT::q1q2_nnlo_r_lz1(z, L)
    +HEFT::q1q2_nnlo_r_lz2(z, L)
    ;
    
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
    
}*/
/*
TEST(LogMufOverMh,N3LO_gg)
{
    const double z=   2.1448782891849743E-002 ;
    const double L=  -2.4079456086518722      ;
    const double falko = HEFT::LEggN3LOregFalko(z, L);
    const double ste = HEFT::LEggNNNLOregSte(z, L);
    cout<<"\nfalko="<<falko;
    cout<<"\nste  ="<<ste<<endl;
    const double err=4e-3;
    EXPECT_LT(abs(falko-ste)/abs(ste),err)
        <<" falko="<<falko
        <<" ste="<<ste
        <<"\t % diff = "<<abs(falko-ste)/abs(ste)
        <<endl;
    
}

TEST(LogMufOverMh,N3LO_gg2)
{
    const double z=   .591448782891849743 ;
    const double L=  2.9324079456086518722      ;
    for (int i=1;i<20;i++)
    {
        double zz=(i)/20.0;
        const double falko = HEFT::LEggN3LOregFalko(zz, L);
        const double ste = HEFT::LEggNNNLOregSte(zz, L);
        cout<<"\nz="<<setw(6)<<zz<<" L="<<setw(8)<<L
        <<"\tfalko="<<setw(9)<<falko
        <<"\tste  ="<<setw(9)<<ste
        <<"\t % diff = "<<setw(12)<<abs(falko-ste)/abs(ste);
        
    }
    
    for (int i=1;i<20;i++)
    {
        double LL=-10.2 + 10.0 * (i)/10.0;
        const double falko = HEFT::LEggN3LOregFalko(z, LL);
        const double ste = HEFT::LEggNNNLOregSte(z, LL);
        cout<<"\nz="<<setw(6)<<fixed<<setprecision(3)<<z<<" L="<<setw(8)<<LL
        <<"\tfalko="<<setw(9)<<falko
        <<"\tste  ="<<setw(9)<<ste
        <<"\t % diff = "<<setw(12)<<scientific<<abs(falko-ste)/abs(ste)<<"\t";
        
    }
    
    const double falko = HEFT::LEggN3LOregFalko(z, L);
    const double ste = HEFT::LEggNNNLOregSte(z, L);
    cout<<"\nfalko="<<falko;
    cout<<"\nste  ="<<ste<<endl;
    const double err=1e-3;
    EXPECT_LT(abs(falko-ste)/abs(ste),err)
    <<" falko="<<falko
    <<" ste="<<ste
    <<"\t % diff = "<<abs(falko-ste)/abs(ste)
    <<endl;
    
}

TEST(LogMufOverMh,N3LO_qg)
{
    const double z=   .591448782891849743 ;
    const double L=  -1.3456      ;
    for (int i=1;i<20;i++)
    {
        double zz=(i)/20.0;
        const double falko = HEFT::LEqgN3LOregFalko(zz, L);
        const double ste = HEFT::LEqgNNNLOregSte(zz, L);
        cout<<"\nz="<<setw(6)<<zz<<" L="<<setw(8)<<L
        <<"\tfalko="<<setw(9)<<falko
        <<"\tste  ="<<setw(9)<<ste
        <<"\t % diff = "<<setw(12)<<abs(falko-ste)/abs(ste);
        
    }
    
    for (int i=1;i<20;i++)
    {
        double LL=-10.+1e-3 + 20.0 * (i)/20.0;
        const double falko = HEFT::LEqgN3LOregFalko(z, LL);
        const double ste = HEFT::LEqgNNNLOregSte(z, LL);
        cout<<"\nz="<<setw(6)<<fixed<<setprecision(3)<<z<<" L="<<setw(8)<<LL
        <<"\tfalko="<<setw(9)<<falko
        <<"\tste  ="<<setw(9)<<ste
        <<"\t % diff = "<<setw(12)<<scientific<<abs(falko-ste)/abs(ste);
        
    }
    
    const double falko = HEFT::LEqgN3LOregFalko(z, L);
    const double ste = HEFT::LEqgNNNLOregSte(z, L);
    cout<<"\nfalko="<<falko;
    cout<<"\nste  ="<<ste<<endl;
    const double err=1e-5;
    EXPECT_LT(abs(falko-ste)/abs(ste),err)
    <<" falko="<<falko
    <<" ste="<<ste
    <<"\t % diff = "<<abs(falko-ste)/abs(ste)
    <<endl;
    
}


TEST(LogMufOverMh,N3LO_qqbar)
{
    const double z=   .591448782891849743 ;
    const double L=  -1.3456      ;
    for (int i=1;i<20;i++)
    {
        double zz=(i)/20.0;
        const double falko = HEFT::LEqqbN3LOregFalko(zz, L);
        const double ste = HEFT::LEqqbNNNLOregSte(zz, L);
        cout<<"\nz="<<setw(6)<<zz<<" L="<<setw(8)<<L
        <<"\tfalko="<<setw(9)<<falko
        <<"\tste  ="<<setw(9)<<ste
        <<"\t % diff = "<<setw(12)<<abs(falko-ste)/abs(ste);
        
    }
    
    for (int i=1;i<20;i++)
    {
        double LL=-10.+1e-3 + 20.0 * (i)/20.0;
        const double falko = HEFT::LEqqbN3LOregFalko(z, LL);
        const double ste = HEFT::LEqqbNNNLOregSte(z, LL);
        cout<<"\nz="<<setw(6)<<fixed<<setprecision(3)<<z<<" L="<<setw(8)<<LL
        <<"\tfalko="<<setw(9)<<falko
        <<"\tste  ="<<setw(9)<<ste
        <<"\t % diff = "<<setw(12)<<scientific<<abs(falko-ste)/abs(ste);
        
    }
    
    const double falko = HEFT::LEqqbN3LOregFalko(z, L);
    const double ste = HEFT::LEqqbNNNLOregSte(z, L);
    cout<<"\nfalko="<<falko;
    cout<<"\nste  ="<<ste<<endl;
    const double err=1e-5;
    EXPECT_LT(abs(falko-ste)/abs(ste),err)
    <<" falko="<<falko
    <<" ste="<<ste
    <<"\t % diff = "<<abs(falko-ste)/abs(ste)
    <<endl;
    
}


TEST(LogMufOverMh,N3LO_qq)
{
    const double z=   .591448782891849743 ;
    const double L=  -1.3456      ;
    for (int i=1;i<20;i++)
    {
        double zz=(i)/20.0;
        const double falko = HEFT::LEqqN3LOregFalko(zz, L);
        const double ste = HEFT::LEqqNNNLOregSte(zz, L);
        cout<<"\nz="<<setw(6)<<zz<<" L="<<setw(8)<<L
        <<"\tfalko="<<setw(9)<<falko
        <<"\tste  ="<<setw(9)<<ste
        <<"\t % diff = "<<setw(12)<<abs(falko-ste)/abs(ste);
        
    }
    
    for (int i=1;i<20;i++)
    {
        double LL=-10.+1e-3 + 20.0 * (i)/20.0;
        const double falko = HEFT::LEqqN3LOregFalko(z, LL);
        const double ste = HEFT::LEqqNNNLOregSte(z, LL);
        cout<<"\nz="<<setw(6)<<fixed<<setprecision(3)<<z<<" L="<<setw(8)<<LL
        <<"\tfalko="<<setw(9)<<falko
        <<"\tste  ="<<setw(9)<<ste
        <<"\t % diff = "<<setw(12)<<scientific<<abs(falko-ste)/abs(ste);
        
    }
    
    const double falko = HEFT::LEqqN3LOregFalko(z, L);
    const double ste = HEFT::LEqqNNNLOregSte(z, L);
    cout<<"\nfalko="<<falko;
    cout<<"\nste  ="<<ste<<endl;
    const double err=1e-5;
    EXPECT_LT(abs(falko-ste)/abs(ste),err)
    <<" falko="<<falko
    <<" ste="<<ste
    <<"\t % diff = "<<abs(falko-ste)/abs(ste)
    <<endl;
    
}


TEST(LogMufOverMh,N3LO_q1q2)
{
    const double z=   .591448782891849743 ;
    const double L=  -1.3456      ;
    for (int i=1;i<20;i++)
    {
        double zz=(i)/20.0;
        const double falko = HEFT::LEq1q2N3LOregFalko(zz, L);
        const double ste = HEFT::LEq1q2NNNLOregSte(zz, L);
        cout<<"\nz="<<setw(6)<<zz<<" L="<<setw(8)<<L
        <<"\tfalko="<<setw(9)<<falko
        <<"\tste  ="<<setw(9)<<ste
        <<"\t % diff = "<<setw(12)<<abs(falko-ste)/abs(ste);
        
    }
    
    for (int i=1;i<20;i++)
    {
        double LL=-10.+1e-3+sqrt(2.)/2. + 20.0 * (i)/20.0;
        const double falko = HEFT::LEq1q2N3LOregFalko(z, LL);
        const double ste = HEFT::LEq1q2NNNLOregSte(z, LL);
        cout<<"\nz="<<setw(6)<<fixed<<setprecision(3)<<z<<" L="<<setw(8)<<LL
        <<"\tfalko="<<setw(9)<<falko
        <<"\tste  ="<<setw(9)<<ste
        <<"\t % diff = "<<setw(12)<<scientific<<abs(falko-ste)/abs(ste);
        
    }
    
    const double falko = HEFT::LEq1q2N3LOregFalko(z, L);
    const double ste = HEFT::LEq1q2NNNLOregSte(z, L);
    cout<<"\nfalko="<<falko;
    cout<<"\nste  ="<<ste<<endl;
    const double err=1e-5;
    EXPECT_LT(abs(falko-ste)/abs(ste),err)
    <<" falko="<<falko
    <<" ste="<<ste
    <<"\t % diff = "<<abs(falko-ste)/abs(ste)
    <<endl;
    
}
*/

/*
// gg channel
// making sure that the series expansions reproduce the full analytic result
// for various values of z: 0.75, 0.25, 0.1
TEST(SeriesExpansionForN3LOgg,coeffOflog5zbar)
{
    int truncation_order = 30;
    const double z=   .591448782891849743 ;
    const double complete = HEFT::gg_n3lo_r_lz5(z,0.0);
    const double series = HEFT::gg_n3lo_r_lz5_series(z,truncation_order);
    
    cout<<"\ncomplete="<<complete;
    cout<<"\nseries (O("<<truncation_order<<"))  ="<<series<<endl;
    const double err=1e-5;
    EXPECT_LT(abs(complete-series)/abs(complete),err)
    <<" complete = "<<complete
    <<" series = "<<series
    <<"\t % diff = "<<abs(complete-series)/abs(complete)
    <<endl;
    
}

TEST(SeriesExpansionForN3LOgg,coeffOflog5zbarCloseToOne)
{
    int truncation_order = 30;
    const double z=   .1591448782891849743 ;
    const double complete = HEFT::gg_n3lo_r_lz5(z,0.0);
    const double series = HEFT::gg_n3lo_r_lz5_series(z,truncation_order);
    
    cout<<"\ncomplete="<<complete;
    cout<<"\nseries (O("<<truncation_order<<"))  ="<<series<<endl;
    const double err=1e-2;
    EXPECT_LT(abs(complete-series)/abs(complete),err)
    <<" complete = "<<complete
    <<" series = "<<series
    <<"\t % diff = "<<abs(complete-series)/abs(complete)
    <<endl;
    
}

TEST(SeriesExpansionForN3LOgg,log5_zbar)
{
    const double z=0.75;
    const double err=1e-10;
    const double res = HEFT::gg_n3lo_r_lz5_series(z,30);
    const double exp = HEFT::gg_n3lo_r_lz5(z,0.0);
    cout<<"\nprecision achieved with 30 terms : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(SeriesExpansionForN3LOgg,log4_zbar)
{
    const double z=0.75;
    const double err=1e-10;
    const double res = HEFT::gg_n3lo_r_lz4_series(z,30);
    const double exp = HEFT::gg_n3lo_r_lz4(z,0.0);
    cout<<"\nprecision achieved with 30 terms : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(SeriesExpansionForN3LOgg,log3_zbar)
{
    const double z=0.75;
    const double err=1e-10;
    const double res = HEFT::gg_n3lo_r_lz3_series(z,30);
    const double exp = HEFT::gg_n3lo_r_lz3(z,0.0);
    cout<<"\nprecision achieved with 30 terms : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(SeriesExpansionForN3LOgg,log5_zbar_z025)
{
    const double z=0.25;
    const double err=1e-3;
    const double res = HEFT::gg_n3lo_r_lz5_series(z,30);
    const double exp = HEFT::gg_n3lo_r_lz5(z,0.0);
    cout<<"\nprecision achieved with 30 terms : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(SeriesExpansionForN3LOgg,log4_zbar_z025)
{
    const double z=0.25;
    const double err=1e-3;
    const double res = HEFT::gg_n3lo_r_lz4_series(z,30);
    const double exp = HEFT::gg_n3lo_r_lz4(z,0.0);
    cout<<"\nprecision achieved with 30 terms : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(SeriesExpansionForN3LOgg,log3_zbar_z025)
{
    const double z=0.25;
    const double err=1e-3;
    const double res = HEFT::gg_n3lo_r_lz3_series(z,30);
    const double exp = HEFT::gg_n3lo_r_lz3(z,0.0);
    cout<<"\nprecision achieved with 30 terms : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(SeriesExpansionForN3LOgg,log5_zbar_z01)
{
    const double z=0.1;
    const double err=1e-1;
    const double res = HEFT::gg_n3lo_r_lz5_series(z,30);
    const double exp = HEFT::gg_n3lo_r_lz5(z,0.0);
    cout<<"\nprecision achieved with 30 terms : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(SeriesExpansionForN3LOgg,log4_zbar_z01)
{
    const double z=0.1;
    const double err=8*1e-1;
    const double res = HEFT::gg_n3lo_r_lz4_series(z,30);
    const double exp = HEFT::gg_n3lo_r_lz4(z,0.0);
    cout<<"\nprecision achieved with 30 terms : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(SeriesExpansionForN3LOgg,log3_zbar_z01)
{
    const double z=0.1;
    const double err=1e-1;
    const double res = HEFT::gg_n3lo_r_lz3_series(z,30);
    const double exp = HEFT::gg_n3lo_r_lz3(z,0.0);
    cout<<"\nprecision achieved with 30 terms : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

*/

// qg channel
// making sure that the series expansions reproduce the full analytic result
// for various values of z: 0.75, 0.25, 0.1 

/*
TEST(SeriesExpansionForN3LOqg,log5_zbar)
{
    const double z=0.75;
    const double err=1e-10;
    const double res = HEFT::qg_n3lo_r_lz5_full_series(z,30);
    const double exp = HEFT::qg_n3lo_r_lz5(z,0.0);
    cout<<"\nprecision achieved with 30 terms : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(SeriesExpansionForN3LOqg,log4_zbar)
{
    const double z=0.75;
    const double err=1e-10;
    const double res = HEFT::qg_n3lo_r_lz4_full_series(z,30);
    const double exp = HEFT::qg_n3lo_r_lz4(z,0.0);
    cout<<"\nprecision achieved with 30 terms : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(SeriesExpansionForN3LOqg,log3_zbar)
{
    const double z=0.75;
    const double err=1e-10;
    const double res = HEFT::qg_n3lo_r_lz3_full_series(z,30);
    const double exp = HEFT::qg_n3lo_r_lz3(z,0.0);
    cout<<"\nprecision achieved with 30 terms : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(SeriesExpansionForN3LOqg,log5_zbar_z025)
{
    const double z=0.25;
    const double err=5*1e-3;
    const double res = HEFT::qg_n3lo_r_lz5_full_series(z,30);
    const double exp = HEFT::qg_n3lo_r_lz5(z,0.0);
    cout<<"\nprecision achieved with 30 terms : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(SeriesExpansionForN3LOqg,log4_zbar_z025)
{
    const double z=0.25;
    const double err=5*1e-3;
    const double res = HEFT::qg_n3lo_r_lz4_full_series(z,30);
    const double exp = HEFT::qg_n3lo_r_lz4(z,0.0);
    cout<<"\nprecision achieved with 30 terms : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(SeriesExpansionForN3LOqg,log3_zbar_z025)
{
    const double z=0.25;
    const double err=5*1e-3;
    const double res = HEFT::qg_n3lo_r_lz3_full_series(z,30);
    const double exp = HEFT::qg_n3lo_r_lz3(z,0.0);
    cout<<"\nprecision achieved with 30 terms : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(SeriesExpansionForN3LOqg,log5_zbar_z01)
{
    const double z=0.1;
    const double err=5*1e-1;
    const double res = HEFT::qg_n3lo_r_lz5_full_series(z,30);
    const double exp = HEFT::qg_n3lo_r_lz5(z,0.0);
    cout<<"\nprecision achieved with 30 terms : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(SeriesExpansionForN3LOqg,log4_zbar_z01)
{
    const double z=0.1;
    const double err=5*1e-1;
    const double res = HEFT::qg_n3lo_r_lz4_full_series(z,30);
    const double exp = HEFT::qg_n3lo_r_lz4(z,0.0);
    cout<<"\nprecision achieved with 30 terms : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(SeriesExpansionForN3LOqg,log3_zbar_z01)
{
    const double z=0.1;
    const double err=5*1e-1;
    const double res = HEFT::qg_n3lo_r_lz3_full_series(z,30);
    const double exp = HEFT::qg_n3lo_r_lz3(z,0.0);
    cout<<"\nprecision achieved with 30 terms : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

*/
/*
TEST(SeriesExpansionForN3LOqqbarAgainstFalko,z03)
{
    const double z=0.3;
    const double err=1e-6;
    double res=0.0;
    for (int i=0;i<6;i++){
        const double curlog = HEFT::qqbar_n3lo_r_lzX_series(i,z,30);
        cout<<fixed<<setprecision(16)<<"\n log(1-z)^"<<i<<" -> "<<curlog<<endl;
        res += curlog;
    }
        
    const double exp = 299.69161198177841;
    cout<<"\nprecision achieved with 30 terms : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}
*/

TEST(LogSeperationForNNLOgg,logzbar_0_logz_0_L0)
{
    const double z=0.1;
    const double err=1e-14;
    const double res = HEFT::nnlo_r_lz0_const(z,0.0);
    const double exp = HEFT::gg_n2lo_lzbar0_lz0_no_log(z);
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(LogSeperationForNNLOgg,logzbar_0_logz_1_L0)
{
    const double z=0.1;
    const double err=1e-15;
    const double res = HEFT::nnlo_r_lz0_logz(z,0.0);
    const double exp = HEFT::gg_n2lo_lzbar0_lz1_no_log(z);
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(LogSeperationForNNLOgg,logzbar_0_logz_2_L0)
{
    const double z=0.1;
    const double err=1e-15;
    const double res = HEFT::nnlo_r_lz0_logz_sq(z,0.0);
    const double exp = HEFT::gg_n2lo_lzbar0_lz2_no_log(z);
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(LogSeperationForNNLOgg,logzbar_0_logz_3_L0)
{
    const double z=0.1;
    const double err=1e-15;
    const double res = HEFT::nnlo_r_lz0_logz_cube(z,0.0);
    const double exp = HEFT::gg_n2lo_lzbar0_lz3_no_log(z);
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}


TEST(LogSeperationForNNLOgg,logzbar_0_Total)
{
    const double z=0.1;
    const double err=1e-12;
    const double L= 0.03;
    const double res = HEFT::nnlo_r_lz0(z,L);
    const double exp = HEFT::gg_n2lo_lzbar0_lz0_no_log(z)
    +HEFT::gg_n2lo_lzbar0_lz1_no_log(z)
    +HEFT::gg_n2lo_lzbar0_lz2_no_log(z)
    +HEFT::gg_n2lo_lzbar0_lz3_no_log(z)
    +HEFT::gg_n2lo_lzbar0_lz0_L_1(z)*L
    +HEFT::gg_n2lo_lzbar0_lz1_L_1(z)*L
    +HEFT::gg_n2lo_lzbar0_lz2_L_1(z)*L
    +HEFT::gg_n2lo_lzbar0_lz0_L_2(z)*L*L
    +HEFT::gg_n2lo_lzbar0_lz1_L_2(z)*L*L
    ;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(LogSeperationForNNLOgg,logzbar_1_Total)
{
    const double z=0.1;
    const double err=1e-12;
    const double L= 45.03;
    const double res = HEFT::nnlo_r_lz1(z,L);
    const double exp = HEFT::gg_n2lo_lzbar1_L_0(z)
    +HEFT::gg_n2lo_lzbar1_L_1(z)*L
    +HEFT::gg_n2lo_lzbar1_L_2(z)*L*L
    ;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(LogSeperationForNNLOgg,logzbar_2_Total)
{
    const double z=0.1;
    const double err=1e-12;
    const double L= 304.03;
    const double res = HEFT::nnlo_r_lz2(z,L);
    const double exp = HEFT::gg_n2lo_lzbar2_L_0(z)
    +HEFT::gg_n2lo_lzbar2_L_1(z)*L
    ;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(LogSeperationForNNLOgg,logzbar_3_Total)
{
    const double z=0.1;
    const double err=1e-12;
    const double L= -30.03;
    const double res = HEFT::nnlo_r_lz3(z,L);
    const double exp = HEFT::gg_n2lo_lzbar3_L_0(z)
    ;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}


TEST(LogSeperationForNNLOgg,Total)
{
    const double z=0.1;
    const double err=1e-12;
    const double L= -30.03;
    const double res = HEFT::nnlo_reg(z,L);
    const double exp = HEFT::gg_n2lo_lzbar0_lz0_no_log(z)
    +HEFT::gg_n2lo_lzbar0_lz1_no_log(z)
    +HEFT::gg_n2lo_lzbar0_lz2_no_log(z)
    +HEFT::gg_n2lo_lzbar0_lz3_no_log(z)
    +HEFT::gg_n2lo_lzbar0_lz0_L_1(z)*L
    +HEFT::gg_n2lo_lzbar0_lz1_L_1(z)*L
    +HEFT::gg_n2lo_lzbar0_lz2_L_1(z)*L
    +HEFT::gg_n2lo_lzbar0_lz0_L_2(z)*L*L
    +HEFT::gg_n2lo_lzbar0_lz1_L_2(z)*L*L
    +HEFT::gg_n2lo_lzbar1_L_0(z)
    +HEFT::gg_n2lo_lzbar1_L_1(z)*L
    +HEFT::gg_n2lo_lzbar1_L_2(z)*L*L
    +HEFT::gg_n2lo_lzbar2_L_0(z)
    +HEFT::gg_n2lo_lzbar2_L_1(z)*L
    +HEFT::gg_n2lo_lzbar3_L_0(z)
    ;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}


TEST(LogSeperationForNNLOgg,Total2)
{
    const double z=0.1;
    const double err=1e-12;
    const double L= 330.03;
    const double res = HEFT::nnlo_reg(z,L);
    const double exp = HEFT::gg_n2lo_lzbar0_lz0_no_log(z)
    +HEFT::gg_n2lo_lzbar0_lz1_no_log(z)
    +HEFT::gg_n2lo_lzbar0_lz2_no_log(z)
    +HEFT::gg_n2lo_lzbar0_lz3_no_log(z)
    +HEFT::gg_n2lo_lzbar0_lz0_L_1(z)*L
    +HEFT::gg_n2lo_lzbar0_lz1_L_1(z)*L
    +HEFT::gg_n2lo_lzbar0_lz2_L_1(z)*L
    +HEFT::gg_n2lo_lzbar0_lz0_L_2(z)*L*L
    +HEFT::gg_n2lo_lzbar0_lz1_L_2(z)*L*L
    +HEFT::gg_n2lo_lzbar1_L_0(z)
    +HEFT::gg_n2lo_lzbar1_L_1(z)*L
    +HEFT::gg_n2lo_lzbar1_L_2(z)*L*L
    +HEFT::gg_n2lo_lzbar2_L_0(z)
    +HEFT::gg_n2lo_lzbar2_L_1(z)*L
    +HEFT::gg_n2lo_lzbar3_L_0(z)
    ;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(LogSeperationForNNLOgg,Total3)
{
    const double z=0.1;
    const double err=1e-12;
    const double L= 0.001;
    const double res = HEFT::nnlo_reg(z,L);
    const double exp = HEFT::gg_n2lo_lzbar0_lz0_no_log(z)
    +HEFT::gg_n2lo_lzbar0_lz1_no_log(z)
    +HEFT::gg_n2lo_lzbar0_lz2_no_log(z)
    +HEFT::gg_n2lo_lzbar0_lz3_no_log(z)
    +HEFT::gg_n2lo_lzbar0_lz0_L_1(z)*L
    +HEFT::gg_n2lo_lzbar0_lz1_L_1(z)*L
    +HEFT::gg_n2lo_lzbar0_lz2_L_1(z)*L
    +HEFT::gg_n2lo_lzbar0_lz0_L_2(z)*L*L
    +HEFT::gg_n2lo_lzbar0_lz1_L_2(z)*L*L
    +HEFT::gg_n2lo_lzbar1_L_0(z)
    +HEFT::gg_n2lo_lzbar1_L_1(z)*L
    +HEFT::gg_n2lo_lzbar1_L_2(z)*L*L
    +HEFT::gg_n2lo_lzbar2_L_0(z)
    +HEFT::gg_n2lo_lzbar2_L_1(z)*L
    +HEFT::gg_n2lo_lzbar3_L_0(z)
    ;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}




TEST(LogSeperationForNNLOqg,logzbar_0)
{
    const double z=0.1;
    const double err=1e-12;
    const double L= -64.32;
    const double res = HEFT::qg_nnlo_r_lz0(z,L);
    const double exp = HEFT::qg_n2lo_r_lz0_L0(z)
    +HEFT::qg_n2lo_r_lz0_L1(z)*L
    +HEFT::qg_n2lo_r_lz0_L2(z)*L*L;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(LogSeperationForNNLOqg,logzbar_1)
{
    const double z=0.1;
    const double err=1e-12;
    const double L= -123.34;
    const double res = HEFT::qg_nnlo_r_lz1(z,L);
    const double exp = HEFT::qg_n2lo_r_lz1_L0(z)
    +HEFT::qg_n2lo_r_lz1_L1(z)*L
    +HEFT::qg_n2lo_r_lz1_L2(z)*L*L;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(LogSeperationForNNLOqg,logzbar_2)
{
    const double z=0.1;
    const double err=1e-12;
    const double L= -45.32;
    const double res = HEFT::qg_nnlo_r_lz2(z,L);
    const double exp = HEFT::qg_n2lo_r_lz2_L0(z)
                      +HEFT::qg_n2lo_r_lz2_L1(z)*L;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(LogSeperationForNNLOqg,logzbar_3)
{
    const double z=0.1;
    const double err=1e-12;
    const double L= -23.45;
    const double res = HEFT::qg_nnlo_r_lz3(z,L);
    const double exp = HEFT::qg_n2lo_r_lz3_L0(z);
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}





TEST(LogSeperationForNNLOqqbar,logzbar_0)
{
    const double z=0.1;
    const double err=1e-12;
    const double L= -12.34;
    const double res = HEFT::qqb_nnlo_r_lz0(z,L);
    const double exp = HEFT::qqb_nnlo_r_lz0_L0(z)
    +HEFT::qqb_nnlo_r_lz0_L1(z)*L
    +HEFT::qqb_nnlo_r_lz0_L2(z)*L*L;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(LogSeperationForNNLOqqbar,logzbar_1)
{
    const double z=0.1;
    const double err=1e-12;
    const double L= -4.54;
    const double res = HEFT::qqb_nnlo_r_lz1(z,L);
    const double exp = HEFT::qqb_nnlo_r_lz1_L0(z)
    +HEFT::qqb_nnlo_r_lz1_L1(z)*L;

    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(LogSeperationForNNLOqqbar,logzbar_2)
{
    const double z=0.1;
    const double err=1e-12;
    const double L= -3.245;
    const double res = HEFT::qqb_nnlo_r_lz2(z,L);
    const double exp = HEFT::qqb_nnlo_r_lz2_L0(z);
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}


TEST(LogSeperationForNNLOqq,logzbar_0)
{
    const double z=0.1;
    const double err=1e-12;
    const double L= 1234.5;
    const double res = HEFT::qq_nnlo_r_lz0(z,L);
    const double exp = HEFT::qq_n2lo_lz0_L0(z)
                      +HEFT::qq_n2lo_lz0_L1(z)*L
                      +HEFT::qq_n2lo_lz0_L2(z)*L*L;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(LogSeperationForNNLOqq,logzbar_1)
{
    const double z=0.1;
    const double err=1e-12;
    const double L= -4.54;
    const double res = HEFT::qq_nnlo_r_lz1(z,L);
    const double exp = HEFT::qq_n2lo_lz1_L0(z)
    +HEFT::qq_n2lo_lz1_L1(z)*L;
    
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(LogSeperationForNNLOqq,logzbar_2)
{
    const double z=0.1;
    const double err=1e-12;
    const double L= -3.245;
    const double res = HEFT::qq_nnlo_r_lz2(z,L);
    const double exp = HEFT::qq_n2lo_lz2_L0(z);
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}









/*
TEST(LogSeperationForNNLOq1q2,logzbar_0)
{
    const double z=0.1;
    const double err=1e-12;
    const double L= 1234.5;
    const double res = HEFT::q1q2_nnlo_r_lz0(z,L);
    const double exp = HEFT::q1q2_n2lo_lz0_L0(z)
    +HEFT::q1q2_n2lo_lz0_L1(z)*L
    +HEFT::q1q2_n2lo_lz0_L2(z)*L*L;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}


TEST(LogSeperationForNNLOq1q2,logzbar_1)
{
    const double z=0.1;
    const double err=1e-12;
    const double L= -4.54;
    const double res = HEFT::q1q2_nnlo_r_lz1(z,L);
    const double exp = HEFT::q1q2_n2lo_lz1_L0(z)
    +HEFT::q1q2_n2lo_lz1_L1(z)*L;
    
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(LogSeperationForNNLOq1q2,logzbar_2)
{
    const double z=0.1;
    const double err=1e-12;
    const double L= -3.245;
    const double res = HEFT::q1q2_nnlo_r_lz2(z,L);
    const double exp = HEFT::q1q2_n2lo_lz2_L0(z);
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}
*/


TEST(LogSeperationForNNLOqq,totalBernhard)
{
    const double z=0.9;
    const double err=1e-12;
    const double L= 0.0;
    const double res = HEFT::qq_n2lo_lz0_L0(z)
    +HEFT::qq_n2lo_lz1_L0(z)
    +HEFT::qq_n2lo_lz2_L0(z)
    +HEFT::qq_n2lo_lz0_L1(z)*L
    +HEFT::qq_n2lo_lz1_L1(z)*L
    +HEFT::qq_n2lo_lz0_L2(z)*L*L;
    
    
    const double exp = 1.5276387979354418;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(LogSeperationForNNLOqq,totalBernhard2)
{
    const double z=0.1;
    const double err=1e-12;
    const double L= 0.0;
    const double res = HEFT::qq_n2lo_lz0_L0(z)
    +HEFT::qq_n2lo_lz1_L0(z)
    +HEFT::qq_n2lo_lz2_L0(z)
    +HEFT::qq_n2lo_lz0_L1(z)*L
    +HEFT::qq_n2lo_lz1_L1(z)*L
    +HEFT::qq_n2lo_lz0_L2(z)*L*L;
    
    
    const double exp = 13.189807539552032;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}
/*
TEST(N3LOggLogDependence,againstBernhard)
{
    const double z=0.9;
    const double err=1e-12;
    const double L= 0.1;
    const double res = HEFT::n3lo_reg_no_Lf(z,37,0.0)
    +HEFT::LEggN3LOregFalko_L1(z)*L
    +HEFT::LEggN3LOregFalko_L2(z)*L*L
    +HEFT::LEggN3LOregFalko_L3(z)*L*L*L;
    
    
    const double exp = -87345.841874745223;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(N3LOggLogDependence,againstBernhardLog1)
{
    const double z=0.2;
    const double err=1e-12;
    const double L= 0.1;
    const double res = HEFT::LEggN3LOregFalko_L1(z);
    
    
    const double exp = -22105.578297257212;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(N3LOggLogDependence,againstBernhardLog2)
{
    const double z=0.2;
    const double err=1e-12;
    const double L= 0.1;
    const double res = HEFT::LEggN3LOregFalko_L2(z);
    
    
    const double exp = -22105.578297257212;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(N3LOggLogDependence,againstBernhardLog3)
{
    for (int i=0;i<20;i++){
        double zz=(i+1)/20.0;
        cout<<"zz="<<zz<<"\t"<<HEFT::LEggN3LOregFalko_L3(zz)<<endl;
    }
    
    
    const double z=0.2;
    const double err=1e-12;
    const double L= 0.1;
    const double res = HEFT::LEggN3LOregFalko_L3(z);
    
    
    const double exp = -22105.578297257212;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(N3LOggreg_new,againstBernhards_new_terms_z_0_9)
{
    for (int i=0;i<19;i++){
        double zz=(i+1)/20.0;
        cout<<"zz="<<zz<<"\t"<<HEFT::n3lo_reg_no_Lf(zz,37,0.0)
        <<"\t"<<HEFT::n3lo_reg_bernhard(zz,37)
        <<"\t"<<(HEFT::n3lo_reg_no_Lf(zz,37,0.0)-HEFT::n3lo_reg_bernhard(zz,37))/HEFT::n3lo_reg_no_Lf(zz,37,0.0)*100.
        <<endl;
    }
    
    
    const double z=0.9;
    const double err=1e-12;
    const double res = HEFT::n3lo_reg_bernhard(z,37);
    
    
    const double exp = HEFT::n3lo_reg_no_Lf(z,37,0.0);
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(N3LOggreg_October15_against_zbar_expansions,Log_zbar_5)
{
    for (int i=0;i<19;i++){
        double zz=(i+1)/20.0;
        double exact=HEFT::gg_n3lo_r_lz5(zz,0);
        double series =HEFT::gg_n3lo_r_lz5_series(zz,37);
        cout<<"zz="<<zz<<"\t"<<series
        <<"\t"<<exact
        <<"\t"<<(exact-series)/exact*100.
        <<endl;
    }
    
    
    const double z=0.9;
    const double err=1e-12;
    const double res = HEFT::gg_n3lo_r_lz5_series(z,37);
    
    
    const double exp = HEFT::gg_n3lo_r_lz5(z,0);
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}


TEST(N3LOggreg_October15_against_zbar_expansions,Log_zbar_4)
{
    for (int i=0;i<19;i++){
        double zz=(i+1)/20.0;
        double exact=HEFT::gg_n3lo_r_lz4(zz,0);
        double series =HEFT::gg_n3lo_r_lz4_series(zz,37);
        cout<<"zz="<<zz<<"\t"<<series
        <<"\t"<<exact
        <<"\t"<<(exact-series)/exact*100.
        <<endl;
    }
    
    
    const double z=0.9;
    const double err=1e-12;
    const double res = HEFT::gg_n3lo_r_lz4_series(z,37);
    
    
    const double exp = HEFT::gg_n3lo_r_lz4(z,0);
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}


TEST(N3LOggreg_October15_against_zbar_expansions,Log_zbar_3)
{
    for (int i=0;i<19;i++){
        double zz=(i+1)/20.0;
        double exact=HEFT::gg_n3lo_r_lz3(zz,0);
        double series =HEFT::gg_n3lo_r_lz3_series(zz,37);
        cout<<"zz="<<zz<<"\t"<<series
        <<"\t"<<exact
        <<"\t"<<(exact-series)/exact*100.
        <<endl;
    }
    
    
    const double z=0.9;
    const double err=1e-12;
    const double res = HEFT::gg_n3lo_r_lz3_series(z,37);
    
    
    const double exp = HEFT::gg_n3lo_r_lz3(z,0);
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}


TEST(N3LOggreg_October15_against_zbar_expansions,Log_zbar_2)
{
    for (int i=0;i<19;i++){
        double zz=(i+1)/20.0;
        double exact=HEFT::gg_n3lo_r_lz2_exact(zz);
        double series =HEFT::gg_n3lo_r_lz2_series(zz,37);
        cout<<"zz="<<zz<<"\t"<<series
        <<"\t"<<exact
        <<"\t"<<(exact-series)/exact*100.
        <<endl;
    }
    
    
    const double z=0.7;
    const double err=1e-7;
    const double res = HEFT::gg_n3lo_r_lz2_series(z,37);
    
    
    const double exp = HEFT::gg_n3lo_r_lz2_exact(z);
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(N3LOggreg_October15_against_zbar_expansions,Log_zbar_1)
{
    for (int i=0;i<19;i++){
        double zz=(i+1)/20.0;
        double exact=HEFT::gg_n3lo_r_lz1_exact(zz);
        double series =HEFT::gg_n3lo_r_lz1_series(zz,37);
        cout<<"zz="<<zz<<"\t"<<series
        <<"\t"<<exact
        <<"\t"<<(exact-series)/exact*100.
        <<endl;
    }
    
    
    const double z=0.7;
    const double err=1e-7;
    const double res = HEFT::gg_n3lo_r_lz1_series(z,37);
    
    
    const double exp = HEFT::gg_n3lo_r_lz1_exact(z);
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}


// qg

TEST(N3LOqgreg_October15_against_zbar_expansions,Log_zbar_2)
{
    for (int i=0;i<19;i++){
        double zz=(i+1)/20.0;
        double exact=HEFT::qg_n3lo_r_lz2_exact(zz);
        double series =HEFT::qg_n3lo_r_lzX_series(2,zz,30);
        cout<<"zz="<<zz<<"\t"<<series
        <<"\t"<<exact
        <<"\t"<<(exact-series)/exact*100.
        <<endl;
    }
    
    
    const double z=0.7;
    const double err=1e-7;
    const double res = HEFT::qg_n3lo_r_lzX_series(2,z,30);
    
    
    const double exp = HEFT::qg_n3lo_r_lz2_exact(z);
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}


TEST(N3LOqgreg_October15_against_zbar_expansions,Log_zbar_1)
{
    for (int i=0;i<19;i++){
        double zz=(i+1)/20.0;
        double exact=HEFT::qg_n3lo_r_lz1_exact(zz);
        double series =HEFT::qg_n3lo_r_lzX_series(1,zz,30);
        cout<<"zz="<<zz<<"\t"<<series
        <<"\t"<<exact
        <<"\t"<<(exact-series)/exact*100.
        <<endl;
    }
    
    
    const double z=0.7;
    const double err=1e-7;
    const double res = HEFT::qg_n3lo_r_lzX_series(1,z,30);
    
    
    const double exp = HEFT::qg_n3lo_r_lz1_exact(z);
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}
*/


/*
TEST(N3LOggreg_October15_me_vs_Bernhard,total_gg)
{
    for (int i=0;i<19;i++){
        double zz=(i+1)/20.0;
        double exact=HEFT::sgg_N3LO_reg(1.-zz);
        double series =HEFT::n3lo_reg_no_Lf(zz,37,0.0);
        cout<<"zz="<<zz<<"\t"<<series
        <<"\t"<<exact
        <<"\t"<<(exact-series)/exact*100.
        <<endl;
    }

    
    const double z=0.7;
    const double err=1e-7;
    const double res = HEFT::n3lo_reg_no_Lf(z,37,0.0);
    
    
    const double exp = HEFT::sgg_N3LO_reg(1.-z);
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}
*/

//TEST(n3Lq1q2simplify,L3)
//{
//    for (int i=0;i<19;i++){
//        double zz=(i+1)/20.0;
//        double exact=HEFT::LEq1q2N3LOregFalko_L3(zz);
//        double new =HEFT::LEq1q2N3LOregFalko_L3new(zz);
//        cout<<"zz="<<zz<<"\t"<<series
//        <<"\t"<<exact
//        <<"\t"<<(exact-series)/exact*100.
//        <<endl;
//    }
//    
//    
//    const double z=0.7;
//    const double err=1e-7;
//    const double res = HEFT::LEq1q2N3LOregFalko_L3(z);
//    
//    
//    const double exp = HEFT::LEq1q2N3LOregFalko_L3new(zz);
//    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
//    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
//}

//: testing zb expansion vs new triple expansion of Bernhard, for reg n3lo gg
//: regions: 1: [0,1/13] , 2: [1/13,0.75], 3: [0.75,1]
class gg_n3lo_reg_triple_expansion_vs_zb_expansion: public ::testing::Test {
    protected:
    virtual void SetUp(){};
    HEFT::N3LORegEvaluator evaluator;
    
};


TEST_F(gg_n3lo_reg_triple_expansion_vs_zb_expansion,DISABLED_Lf_zero_region_1) {
    const double err=1e-7;
    const double z=0.03;
    double exp = HEFT::n3lo_reg_no_Lf(z,37,0.0);
    double res = evaluator.n3lo_reg_complete_gg(z,0);
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST_F(gg_n3lo_reg_triple_expansion_vs_zb_expansion,Lf_zero_region_2) {
    const double err=1e-7;
    const double z=0.43;
    double exp = HEFT::n3lo_reg_no_Lf(z,37,0.0);
    double res = evaluator.n3lo_reg_complete_gg(z,0);
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST_F(gg_n3lo_reg_triple_expansion_vs_zb_expansion,Lf_zero_region_3) {
    const double err=1e-7;
    const double z=0.93;
    double exp = HEFT::n3lo_reg_no_Lf(z,37,0.0);
    double res = evaluator.n3lo_reg_complete_gg(z,0);
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST_F(gg_n3lo_reg_triple_expansion_vs_zb_expansion,DISABLED_Lf_zero_region_all) {
    const double err=1e-7;
    double total_loss_prec=0.0;
    for (int i=1;i<1000;i++) {
        const double z=0.001 * double(i);
        double exp = HEFT::n3lo_reg_no_Lf(z,37,0.0);
        double res = evaluator.n3lo_reg_complete_gg(z,0);
        total_loss_prec += abs(res-exp)/abs(exp);
        cout<<setprecision(16)<<"\n{" << z << "," << exp << "," << res << "},";
    }
    EXPECT_LT(total_loss_prec,err)<<" res="<<total_loss_prec<<" expected = "<<err<<endl;
}



TEST_F(gg_n3lo_reg_triple_expansion_vs_zb_expansion,Lf_1_region_all) {
    const double err=1e-7;
    double total_loss_prec=0.0;
    for (int i=0;i<100;i++) {
        const double z=0.1 * double(i)/10.0+0.009;
        double exp = HEFT::LEggN3LOregFalko_L1(z);
        double res = -evaluator.n3lo_reg_complete_gg(z,1);
        total_loss_prec += abs(res-exp)/abs(exp);
        cout<<"\n z = " << z << " old = " << exp << " new = " << res << " precision achieved : "<<abs(res-exp)/abs(exp)<<endl;
        }
    EXPECT_LT(total_loss_prec,err)<<" res="<<total_loss_prec<<" expected = "<<err<<endl;
}

TEST_F(gg_n3lo_reg_triple_expansion_vs_zb_expansion,Lf_2_region_all) {
    const double err=1e-7;
    double total_loss_prec=0.0;
    for (int i=0;i<100;i++) {
        const double z=0.1 * double(i)/10.0+0.009;
        double exp = HEFT::LEggN3LOregFalko_L2(z);
        double res = evaluator.n3lo_reg_complete_gg(z,2);
        total_loss_prec += abs(res-exp)/abs(exp);
        cout<<"\n z = " << z << " old = " << exp << " new = " << res << " precision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    }
    EXPECT_LT(total_loss_prec,err)<<" res="<<total_loss_prec<<" expected = "<<err<<endl;
}

TEST_F(gg_n3lo_reg_triple_expansion_vs_zb_expansion,Lf_3_region_all) {
    const double err=1e-7;
    double total_loss_prec=0.0;
    for (int i=0;i<100;i++) {
        const double z=0.1 * double(i)/10.0+0.009;
        double exp = HEFT::LEggN3LOregFalko_L3(z);
        double res = -evaluator.n3lo_reg_complete_gg(z,3);
        total_loss_prec += abs(res-exp)/abs(exp);
        cout<<"\n z = " << z << " old = " << exp << " new = " << res << " precision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    }
    EXPECT_LT(total_loss_prec,err)<<" res="<<total_loss_prec<<" expected = "<<err<<endl;
}

class qg_n3lo_reg_triple_expansion_vs_zb_expansion: public ::testing::Test {
    protected:
    virtual void SetUp(){};
    HEFT::N3LORegEvaluator evaluator;
    
};

TEST_F(qg_n3lo_reg_triple_expansion_vs_zb_expansion,DISABLED_Lf_zero_region_all) {
    const double err=1e-7;
    double total_loss_prec=0.0;
    for (int i=1;i<100;i++) {
        const double z=0.01 * double(i);
        double exp = HEFT::qg_n3lo_reg_no_Lf(z,37,0.0);
        double res = evaluator.n3lo_reg_complete_qg(z,0);
        total_loss_prec += abs(res-exp)/abs(exp);
        cout<<setprecision(16)<<"\n{" << z << "," << exp << "," << res << "},";
    }
    EXPECT_LT(total_loss_prec,err)<<" res="<<total_loss_prec<<" expected = "<<err<<endl;
}

class qqbar_n3lo_reg_triple_expansion_vs_zb_expansion: public ::testing::Test {
    protected:
    virtual void SetUp(){};
    HEFT::N3LORegEvaluator evaluator;
    
};

TEST_F(qqbar_n3lo_reg_triple_expansion_vs_zb_expansion,DISABLED_Lf_zero_region_all) {
    const double err=1e-7;
    double total_loss_prec=0.0;
    for (int i=1;i<100;i++) {
        const double z=0.01 * double(i);
        double exp = HEFT::qqb_n3lo_reg_no_Lf(z,37);
        double res = evaluator.n3lo_reg_complete_qqbar(z,0);
        total_loss_prec += abs(res-exp)/abs(exp);
        cout<<setprecision(16)<<"\n{" << z << "," << exp << "," << res << "},";
    }
    EXPECT_LT(total_loss_prec,err)<<" res="<<total_loss_prec<<" expected = "<<err<<endl;
}

class qq_n3lo_reg_triple_expansion_vs_zb_expansion: public ::testing::Test {
    protected:
    virtual void SetUp(){};
    HEFT::N3LORegEvaluator evaluator;
    
};

TEST_F(qq_n3lo_reg_triple_expansion_vs_zb_expansion,DISABLED_Lf_zero_region_all) {
    const double err=1e-7;
    double total_loss_prec=0.0;
    for (int i=1;i<100;i++) {
        const double z=0.01 * double(i);
        double exp = HEFT::qq_n3lo_reg_no_Lf(z,37);
        double res = evaluator.n3lo_reg_complete_qq(z,0);
        total_loss_prec += abs(res-exp)/abs(exp);
        cout<<setprecision(16)<<"\n{" << z << "," << exp << "," << res << "},";
    }
    EXPECT_LT(total_loss_prec,err)<<" res="<<total_loss_prec<<" expected = "<<err<<endl;
}

class q1q2_n3lo_reg_triple_expansion_vs_zb_expansion: public ::testing::Test {
    protected:
    virtual void SetUp(){};
    HEFT::N3LORegEvaluator evaluator;
    
};

TEST_F(q1q2_n3lo_reg_triple_expansion_vs_zb_expansion,DISABLED_Lf_zero_region_all) {
    const double err=1e-7;
    double total_loss_prec=0.0;
    for (int i=1;i<100;i++) {
        const double z=0.01 * double(i);
        double exp = HEFT::q1q2_n3lo_reg_no_Lf(z,37);
        double res = evaluator.n3lo_reg_complete_q1q2(z,0);
        total_loss_prec += abs(res-exp)/abs(exp);
        cout<<setprecision(16)<<"\n{" << z << "," << exp << "," << res << "},";
    }
    EXPECT_LT(total_loss_prec,err)<<" res="<<total_loss_prec<<" expected = "<<err<<endl;
}



TEST(N3LO_delta,Lf_0)
{
   const double err=1e-7;
    const double res = HEFT::n_N3LO_delta();
    const double exp = 1124.308887494467;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(N3LO_delta,Lf_1)
{
    const double err=1e-7;
    const double res = HEFT::n_N3LO_delta_L();
    const double exp = -316.3371547476737;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(N3LO_delta,Lf_2)
{
    const double err=1e-7;
    const double res = HEFT::n_N3LO_delta_L2();
    const double exp = -746.7559472635107;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(N3LO_delta,Lf_3)
{
    const double err=1e-7;
    const double res = HEFT::n_N3LO_delta_L3();
    const double exp = -143.2983223337546;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}


TEST(N3LO_D0,Lf_0)
{
    const double err=1e-7;
    const double res = HEFT::n_N3LO_D0();
    const double exp = 1466.478272427120;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(N3LO_D0,Lf_1)
{
    const double err=1e-7;
    const double res = HEFT::n_N3LO_D0_L();
    const double exp = 3031.043368878362;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}
TEST(N3LO_D0,Lf_2)
{
    const double err=1e-7;
    const double res = HEFT::n_N3LO_D0_L2();
    const double exp = 1253.083450946629;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}
TEST(N3LO_D0,Lf_3)
{
    const double err=1e-7;
    const double res = HEFT::n_N3LO_D0_L3();
    const double exp = 170.3056569973862;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}


TEST(N3LO_D1,Lf_0)
{
    const double err=1e-7;
    const double res = HEFT::n_N3LO_D1();
    const double exp = -6062.086737756724;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(N3LO_D1,Lf_1)
{
    const double err=1e-7;
    const double res = HEFT::n_N3LO_D1_L();
    const double exp = -6732.891296207376;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}
TEST(N3LO_D1,Lf_2)
{
    const double err=1e-7;
    const double res = HEFT::n_N3LO_D1_L2();
    const double exp = -1252.792579643143;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}
TEST(N3LO_D1,Lf_3)
{
    const double err=1e-7;
    const double res = HEFT::n_N3LO_D1_L3();
    const double exp = 69.00000000000000;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(N3LO_D2,Lf_0)
{
    const double err=1e-7;
    const double res = HEFT::n_N3LO_D2();
    const double exp = 7116.015301709462;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(N3LO_D2,Lf_1)
{
    const double err=1e-7;
    const double res = HEFT::n_N3LO_D2_L();
    const double exp = 2736.543796945111;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}
TEST(N3LO_D2,Lf_2)
{
    const double err=1e-7;
    const double res = HEFT::n_N3LO_D2_L2();
    const double exp = -310.5000000000000;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}
TEST(N3LO_D2,Lf_3)
{
    const double err=1e-7;
    const double res = HEFT::n_N3LO_D2_L3();
    const double exp = -108.0000000000000;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(N3LO_D3,Lf_0)
{
    const double err=1e-7;
    const double res = HEFT::n_N3LO_D3();
    const double exp = -1824.362531296741;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(N3LO_D3,Lf_1)
{
    const double err=1e-7;
    const double res = HEFT::n_N3LO_D3_L();
    const double exp = 460.0000000000000;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}
TEST(N3LO_D3,Lf_2)
{
    const double err=1e-7;
    const double res = HEFT::n_N3LO_D3_L2();
    const double exp = 432.0000000000000;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(N3LO_D4,Lf_0)
{
    const double err=1e-7;
    const double res = HEFT::n_N3LO_D4();
    const double exp = -230.0000000000000;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(N3LO_D4,Lf_1)
{
    const double err=1e-7;
    const double res = HEFT::n_N3LO_D4_L();
    const double exp = -540.0000000000000;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}

TEST(N3LO_D5,Lf_0)
{
    const double err=1e-7;
    const double res = HEFT::n_N3LO_D5();
    const double exp = 216.0000000000000;
    cout<<"\nprecision achieved : "<<abs(res-exp)/abs(exp)<<endl;
    EXPECT_LT(abs((res-exp)/exp),err)<<" res="<<res<<" expected = "<<exp<<endl;
}
int main(int argc, char**argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return  RUN_ALL_TESTS();
    return 0;
}




















