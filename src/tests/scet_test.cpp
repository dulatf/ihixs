/** testing UserInterface.*:
 *
 * Achilleas Lazopoulos, lazopoli@phys.ethz.ch
 */

#include <iostream>
#include <cmath>



using namespace std;

#include "gtest/gtest.h"
#include "scet_resum.h"

class TestAgainstBernhardsCode: public ::testing::Test {
protected:
    virtual void SetUp()
    {
        UI.m_higgs = 125.5;
        
        UI.top_scheme = "on-shell";
        UI.mt_on_shell = 173.4 ;//# 173.4 pdg
        
        UI.bottom_scheme="on-shell";
        UI.mb_on_shell = 4.75;// #  from pdg http://pdg.lbl.gov/2014/listings/rpp2014-list-b-quark.pdf
        
        UI.charm_scheme = "on-shell";
        UI.mc_on_shell = 1.40 ;//# 1.67 from pdg http://pdg.lbl.gov/2014/listings/rpp2014-list-c-quark.pdf
        
        
        
        
        UI.mur=125.5;
        UI.muf=125.5;
        UI.Etot = 13000.;
        UI.pdf_set = "NNPDF30_nnlo_as_0118";
        UI.qcd_perturbative_order = "N3LO";
        
    }
    
    UserInterface UI;
    InputParameters IP;
};


TEST_F(TestAgainstBernhardsCode,Order2AtMh)
{
    UI.mur=UI.m_higgs;
    UI.muf=UI.m_higgs;
    
    IP.Configure(UI);
    
    double res = SCET::ScetResummation(UI,&IP,2);
    double exp = 5.53676;
    double err=1e-4;
    
    
    
    EXPECT_LT(abs(res-exp)/abs(exp),err)<<"res="<<res<<"| expecting="<<exp<<"|"<<endl;
}

TEST_F(TestAgainstBernhardsCode,Order2AtMhOver2)
{
    UI.mur=UI.m_higgs/2.;
    UI.muf=UI.m_higgs/2.;
    
    IP.Configure(UI);
    
    double res = SCET::ScetResummation(UI,&IP,2);
    double exp = 0.79657;
    double err=1e-4;
    
    
    
    EXPECT_LT(abs(res-exp)/abs(exp),err)<<"res="<<res<<"| expecting="<<exp<<"|"<<endl;
}

TEST_F(TestAgainstBernhardsCode,Order2AtMhOver4)
{
    UI.mur=UI.m_higgs/4.;
    UI.muf=UI.m_higgs/4.;
    
    IP.Configure(UI);
    
    double res = SCET::ScetResummation(UI,&IP,2);
    double exp = -0.90432;
    double err=1e-4;
    
    
    
    EXPECT_LT(abs(res-exp)/abs(exp),err)<<"res="<<res<<"| expecting="<<exp<<"|"<<endl;
}

TEST_F(TestAgainstBernhardsCode,Order2AtMhTimes2)
{
    UI.mur=UI.m_higgs*2.;
    UI.muf=UI.m_higgs*2.;
    
    IP.Configure(UI);
    
    double res = SCET::ScetResummation(UI,&IP,2);
    double exp = 12.59946;
    double err=1e-4;
    
    
    
    EXPECT_LT(abs(res-exp)/abs(exp),err)<<"res="<<res<<"| expecting="<<exp<<"|"<<endl;
}


TEST_F(TestAgainstBernhardsCode,Order3AtMh)
{
    UI.mur=UI.m_higgs;
    UI.muf=UI.m_higgs;
    
    IP.Configure(UI);
    
    double res = SCET::ScetResummation(UI,&IP,3);
    double exp = 1.53590;
    double err=1e-4;
    
    
    
    EXPECT_LT(abs(res-exp)/abs(exp),err)<<"res="<<res<<"| expecting="<<exp<<"|"<<endl;
}

/*
TEST_F(TestAgainstBernhardsCode,Order3AtMhOver2)
{
    UI.mur=UI.m_higgs/2.;
    UI.muf=UI.m_higgs/2.;
    
    IP.Configure(UI);
    
    double res = SCET::ScetResummation(UI,&IP,3);
    double exp =  -0.00572;
    double err=1e-3;
    
    
    
    EXPECT_LT(abs(res-exp)/abs(exp),err)<<"res="<<res<<"| expecting="<<exp<<"|"<<endl;
}
 */
/*
TEST_F(TestAgainstBernhardsCode,Order3AtMhOver4)
{
    UI.mur=UI.m_higgs/4.;
    UI.muf=UI.m_higgs/4.;
    
    IP.Configure(UI);
    
    double res = SCET::ScetResummation(UI,&IP,3);
    double exp = 0.22767;
    double err=1e-4;
    
    
    
    EXPECT_LT(abs(res-exp)/abs(exp),err)<<"res="<<res<<"| expecting="<<exp<<"|"<<endl;
}

TEST_F(TestAgainstBernhardsCode,Order3AtMhTimes2)
{
    UI.mur=UI.m_higgs*2.;
    UI.muf=UI.m_higgs*2.;
    
    IP.Configure(UI);
    
    double res = SCET::ScetResummation(UI,&IP,3);
    double exp = 5.09697;
    double err=1e-4;
    
    
    
    EXPECT_LT(abs(res-exp)/abs(exp),err)<<"res="<<res<<"| expecting="<<exp<<"|"<<endl;
}

TEST_F(TestAgainstBernhardsCode,MurDiffThanMuf)
{
    UI.mur=UI.m_higgs*2.;
    UI.muf=UI.m_higgs/2.;
    
    IP.Configure(UI);
    
    double res = SCET::ScetResummation(UI,&IP,3);
    double exp = 0.0;
    double err=1e-4;
    
    
    
    EXPECT_EQ(res,exp)<<"res="<<res<<"| expecting="<<exp<<"|"<<endl;
}*/

int main(int argc, char**argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return  RUN_ALL_TESTS();
    return 0;
}




















