/** testing UserInterface.*:
 *
 * Achilleas Lazopoulos, lazopoli@phys.ethz.ch
 */

#include <iostream>
#include <cmath>



using namespace std;

#include "gtest/gtest.h"
#include "mt_expansion_qg.h"
#include "mt_expansion_gg.h"
#include "inclusive_higgs_eft.h"
#include "inclusive_mt_expansion.h"
#include "higgs_eft.h"
#include "math.h"
#include "wilson_coefficients.h"


TEST(HighEnergyLimitTest,linear_interpolator_on_grid_point)
{
    double y[23]={1, 4, 9, 16, 25, 36, 49, 64, 81, 100, 121, 144, 169, 196, 225, 256, \
        289, 324, 361, 400, 441, 484, 529};
    const double res = MTEXP::linear_interpolate(y,6.0);
    const double exp = 121.;
    const double err = 1e-15;
    EXPECT_LT(abs((res-exp)/exp),err);
}

TEST(HighEnergyLimitTest,linear_interpolator_in_between)
{
    double y[23]={1, 4, 9, 16, 25, 36, 49, 64, 81, 100, 121, 144, 169, 196, 225, 256, \
        289, 324, 361, 400, 441, 484, 529};
    const double res = MTEXP::linear_interpolate(y,6.3);
    const double exp = 134.8;
    const double err = 1e-15;
    EXPECT_LT(abs((res-exp)/exp),err);
}


class TestAgainstEFTLargeMt: public ::testing::Test {
protected:
    virtual void SetUp()
    {
        UI.SetOption("m_higgs" , 125.);
        
        UI.SetOption("top_scheme" , "on-shell");
        UI.SetOption("mt_on_shell",17000.9) ;//# 173.4 pdg
        
        UI.SetOption("bottom_scheme" , "on-shell");
        UI.SetOption("mb_on_shell" , 4.75);// #  from pdg http://pdg.lbl.gov/2014/listings/rpp2014-list-b-quark.pdf
        
        UI.SetOption("charm_scheme" , "on-shell");
        UI.SetOption("mc_on_shell" , 1.40) ;//# 1.67 from pdg http://pdg.lbl.gov/2014/listings/rpp2014-list-c-quark.pdf
        
        
        
        
        UI.SetOption("mur",25.);
        UI.SetOption("muf",25.);
        UI.SetOption("Etot" , 14000.);
        
        
        UI.SetOption("pdf_set" , "PDF4LHC15_nnlo_100");
        UI.SetOption("qcd_perturbative_order" , "N3LO");
        
        
        
    }
    
    UserInterface UI;
};


 TEST_F(TestAgainstEFTLargeMt,ggdelta)
 {
     InclusiveHiggsEFT* eft = new InclusiveHiggsEFT(UI);
     eft->SetUpContributions();
     eft->EvaluateSingleTerm("gg_delta_lo");
     ResultPair eftres=eft->RescalingCoefficient()*eft->GiveChannel("gg")->TermPtr("gg_delta_lo")->QCDResult().term_of_order(4);
     
     InclusiveHiggsMtExpansion* mt_expansion = new InclusiveHiggsMtExpansion(UI);
     mt_expansion->SetUpContributions();
     
     mt_expansion->EvaluateSingleTerm("gg_delta_nnlo_mt_exp");
     ResultPair mtres=mt_expansion->GiveChannel("gg")->TermPtr("gg_delta_nnlo_mt_exp")->QCDResult().term_of_order(4);
     
 
     const double err=1e-3;
     EXPECT_LT(abs(eftres.val()-mtres.val())/abs(eftres.val()),err)<<"eft="<<eftres.val()<<"| mtexp="<<mtres.val()<<"|"<<endl;
 }

TEST_F(TestAgainstEFTLargeMt,D0)
{
    InclusiveHiggsEFT* eft = new InclusiveHiggsEFT(UI);
    eft->SetUpContributions();
    eft->EvaluateSingleTerm("gg_D0");
    ResultPair eftres=eft->RescalingCoefficient()*eft->GiveChannel("gg")->TermPtr("gg_D0")->QCDResult().term_of_order(4);
    
    InclusiveHiggsMtExpansion* mt_expansion = new InclusiveHiggsMtExpansion(UI);
    mt_expansion->SetUpContributions();
    
    mt_expansion->EvaluateSingleTerm("gg_D0_nnlo_mt_exp");
    ResultPair mtres=mt_expansion->GiveChannel("gg")->TermPtr("gg_D0_nnlo_mt_exp")->QCDResult().term_of_order(4);
    
    
    const double err=1e-3;
    EXPECT_LT(abs(eftres.val()-mtres.val())/abs(eftres.val()),err)<<"eft="<<eftres.val()<<"| mtexp="<<mtres.val()<<"|"<<endl;
}

TEST_F(TestAgainstEFTLargeMt,D1)
{
    InclusiveHiggsEFT* eft = new InclusiveHiggsEFT(UI);
    eft->SetUpContributions();
    eft->EvaluateSingleTerm("gg_D1");
    ResultPair eftres=eft->RescalingCoefficient()*eft->GiveChannel("gg")->TermPtr("gg_D1")->QCDResult().term_of_order(4);
    
    InclusiveHiggsMtExpansion* mt_expansion = new InclusiveHiggsMtExpansion(UI);
    mt_expansion->SetUpContributions();
    
    mt_expansion->EvaluateSingleTerm("gg_D1_nnlo_mt_exp");
    ResultPair mtres=mt_expansion->GiveChannel("gg")->TermPtr("gg_D1_nnlo_mt_exp")->QCDResult().term_of_order(4);
    
    
    const double err=1e-3;
    EXPECT_LT(abs(eftres.val()-mtres.val())/abs(eftres.val()),err)<<"eft="<<eftres.val()<<"| mtexp="<<mtres.val()<<"|"<<endl;
}


TEST_F(TestAgainstEFTLargeMt,D2)
{
    InclusiveHiggsEFT* eft = new InclusiveHiggsEFT(UI);
    eft->SetUpContributions();
    eft->EvaluateSingleTerm("gg_D2");
    ResultPair eftres=eft->RescalingCoefficient()*eft->GiveChannel("gg")->TermPtr("gg_D2")->QCDResult().term_of_order(4);
    
    InclusiveHiggsMtExpansion* mt_expansion = new InclusiveHiggsMtExpansion(UI);
    mt_expansion->SetUpContributions();
    
    mt_expansion->EvaluateSingleTerm("gg_D2_nnlo_mt_exp");
    ResultPair mtres=mt_expansion->GiveChannel("gg")->TermPtr("gg_D2_nnlo_mt_exp")->QCDResult().term_of_order(4);
    
    
    
    const double err=1e-3;
    EXPECT_LT(abs(eftres.val()-mtres.val())/abs(eftres.val()),err)<<"eft="<<eftres.val()<<"| mtexp="<<mtres.val()<<"|"<<endl;
}

TEST_F(TestAgainstEFTLargeMt,D3)
{
    InclusiveHiggsEFT* eft = new InclusiveHiggsEFT(UI);
    eft->SetUpContributions();
    eft->EvaluateSingleTerm("gg_D3");
    ResultPair eftres=eft->RescalingCoefficient()*eft->GiveChannel("gg")->TermPtr("gg_D3")->QCDResult().term_of_order(4);
    
    InclusiveHiggsMtExpansion* mt_expansion = new InclusiveHiggsMtExpansion(UI);
    mt_expansion->SetUpContributions();
    
    mt_expansion->EvaluateSingleTerm("gg_D3_nnlo_mt_exp");
    ResultPair mtres=mt_expansion->GiveChannel("gg")->TermPtr("gg_D3_nnlo_mt_exp")->QCDResult().term_of_order(4);
  
    
    
    const double err=1e-3;
    EXPECT_LT(abs(eftres.val()-mtres.val())/abs(eftres.val()),err)<<"eft="<<eftres.val()<<"| mtexp="<<mtres.val()<<"|"<<endl;
}

TEST_F(TestAgainstEFTLargeMt,DISABLED_reg)
{
    UI.SetOption("mt_on_shell" , 216.0) ;//# 173.4 pdg
    InclusiveHiggsEFT* eft = new InclusiveHiggsEFT(UI);
    eft->SetUpContributions();
    eft->EvaluateSingleTerm("gg_reg_n2lo");
    eft->EvaluateSingleTerm("gg_reg_nlo");
    ResultPair eftres=eft->RescalingCoefficient() * (
                eft->GiveChannel("gg")->TermPtr("gg_reg_n2lo")->QCDResult().term_of_order(4)
                +eft->GiveChannel("gg")->TermPtr("gg_reg_nlo")->QCDResult().term_of_order(4));
    
    InclusiveHiggsMtExpansion* mt_expansion = new InclusiveHiggsMtExpansion(UI);
    mt_expansion->SetUpContributions();
    
    mt_expansion->EvaluateSingleTerm("gg_reg_nnlo_mt_exp");
    ResultPair mtres=mt_expansion->GiveChannel("gg")->TermPtr("gg_reg_nnlo_mt_exp")->QCDResult().term_of_order(4);
    
    
    
    
    const double err=1e-3;
    EXPECT_LT(abs(eftres.val()-mtres.val())/abs(eftres.val()),err)<<"eft="<<eftres.val()<<"| mtexp="<<mtres.val()<<"|"<<endl
    <<"***: It is reasonable that this test fails: it effectively compares the complete EFT with the z-bar expanded EFT"<<endl;
}

TEST_F(TestAgainstEFTLargeMt,DISABLED_qg_reg)
{
    UI.SetOption("mt_on_shell" , 170.9) ;//# 173.4 pdg
    InclusiveHiggsEFT* eft = new InclusiveHiggsEFT(UI);
    eft->SetUpContributions();
    eft->EvaluateSingleTerm("qg_reg_n2lo");
    eft->EvaluateSingleTerm("qg_reg_nlo");
    ResultPair eftres=eft->RescalingCoefficient()* (
                eft->GiveChannel("qg")->TermPtr("qg_reg_n2lo")->QCDResult().term_of_order(4)
               +eft->GiveChannel("qg")->TermPtr("qg_reg_nlo")->QCDResult().term_of_order(4));
    
    InclusiveHiggsMtExpansion* mt_expansion = new InclusiveHiggsMtExpansion(UI);
    mt_expansion->SetUpContributions();
    
    mt_expansion->EvaluateSingleTerm("qg_reg_nnlo_mt_exp");
    ResultPair mtres=mt_expansion->GiveChannel("qg")->TermPtr("qg_reg_nnlo_mt_exp")->QCDResult().term_of_order(4);
    
    
    
    
    const double err=1e-3;
    EXPECT_LT(abs(eftres.val()-mtres.val())/abs(eftres.val()),err)<<"eft="<<eftres.val()<<"| mtexp="<<mtres.val()<<"|"<<endl
    <<"***: It is reasonable that this test fails: it effectively compares the complete EFT with the z-bar expanded EFT"<<endl;
}


TEST(QuarkGluonSingleZValue,z0p1)
{
    const double mh=125.;
    const double mt=175.;
    const double z=0.1;
    const double rho = pow(mh/mt,2.);
    HiggsMtExpansion_qg_nnlo_reg  mtqgreg(rho);
    
    double res=mtqgreg.z_times_reg_L0(z);
    
    
    const double exp = 19.294;
    const double err = 1e-4;
    EXPECT_LT(abs((res-exp)/exp),err)<<"res = "<<res<<"\t exp = "<<exp<<endl;
}

TEST(QuarkGluonSingleZValueAainstEFT,z0p345)
{
    //const double mh=125.;
    const double z=0.345;
    const double rho = 0.0001;//pow(mh/mt,2.);
    HiggsMtExpansion_qg_nnlo_reg mtqgreg(rho);
    
    double res=mtqgreg.z_times_reg_L0(z)/z;
    
    
    const double exp = HEFT::qg_nnlo_reg(z,0.0)+11./2.*HEFT::qg_nlo_reg(z,0.0);
    const double err = 1e-4;
    EXPECT_LT(abs((res-exp)/exp),err)<<"res = "<<res<<"\t exp = "<<exp<<endl;
}


TEST(QuarkGluonSingleZValueAainstEFT,zscan)
{
    //const double mh=125.;
    const double z=0.345;
    const double rho = 0.0001;//pow(mh/mt,2.);
    HiggsMtExpansion_qg_nnlo_reg mtqgreg(rho);
    
    for (int i=0;i<18;i++){
        double zz=0.05+i*0.05;
        cout<<"\nzz="<<zz<<"\tmt="<<mtqgreg.z_times_reg_L0(zz)/zz
        <<"\teft="<<HEFT::qg_nnlo_reg(zz,0.0)+11./2.*HEFT::qg_nlo_reg(zz,0.0)<<endl;
    }
    
    double res=mtqgreg.z_times_reg_L0(z)/z;
    
    
    const double exp = HEFT::qg_nnlo_reg(z,0.0)+11./2.*HEFT::qg_nlo_reg(z,0.0);
    const double err = 1e-4;
    EXPECT_LT(abs((res-exp)/exp),err)<<"res = "<<res<<"\t exp = "<<exp<<endl
    <<"***: It is reasonable that this test fails: it effectively compares the complete EFT with the z-bar expanded EFT"<<endl;
}

int main(int argc, char**argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return  RUN_ALL_TESTS();
    return 0;
}


class TestAgainstMathematica: public ::testing::Test {
protected:
    virtual void SetUp()
    {
        UI.SetOption("m_higgs" , 125.);
        
        UI.SetOption("top_scheme" , "on-shell");
        UI.SetOption("mt_on_shell",170.) ;//# 173.4 pdg
        
        UI.SetOption("bottom_scheme" , "on-shell");
        UI.SetOption("mb_on_shell" , 4.75);// #  from pdg http://pdg.lbl.gov/2014/listings/rpp2014-list-b-quark.pdf
        
        UI.SetOption("charm_scheme" , "on-shell");
        UI.SetOption("mc_on_shell" , 1.40) ;//# 1.67 from pdg http://pdg.lbl.gov/2014/listings/rpp2014-list-c-quark.pdf
        
        
        
        
        UI.SetOption("mur",25.);
        UI.SetOption("muf",25.);
        UI.SetOption("Etot" , 13000.);
        
        
        UI.SetOption("pdf_set" , "NNPDF30_nnlo_as_0118");
        UI.SetOption("qcd_perturbative_order" , "N3LO");
        
        rho = pow(UI.giveDouble("m_higgs")/UI.giveDouble("mt_on_shell"),2.);
        lmh = 2.*log(UI.giveDouble("mur")/UI.giveDouble("m_higgs"));
                     lmtop=2.*log(UI.giveDouble("mur")/UI.giveDouble("mt_on_shell"));

    
       
    }
    
    UserInterface UI;
    double rho,lmh,lmtop;
};

TEST_F(TestAgainstMathematica,quark_gluon_reg_L0)
{
    const double z=0.345;
    HiggsMtExpansion_qg_nnlo_reg mtqgreg(rho);
    
    double res=mtqgreg.z_times_reg_L0(z)/z;
    
    
    const double exp = -17.15443252946232;
    const double err = 1e-4;
    EXPECT_LT(abs((res-exp)/exp),err)<<"res = "<<res<<"\t exp = "<<exp<<endl;
}

TEST_F(TestAgainstMathematica,quark_gluon_reg_L1)
{
    const double z=0.345;
    HiggsMtExpansion_qg_nnlo_reg mtqgreg(rho);
    
    double res=mtqgreg.z_times_reg_L1(z)/z;
    
    
    const double exp = -76.41931650942522;
    const double err = 1e-4;
    EXPECT_LT(abs((res-exp)/exp),err)<<"res = "<<res<<"\t exp = "<<exp<<endl;
}




TEST_F(TestAgainstMathematica,quark_gluon_reg_L2)
{
    const double z=0.345;
    HiggsMtExpansion_qg_nnlo_reg mtqgreg(rho);
    
    double res=mtqgreg.z_times_reg_L2(z)/z;
    
    
    const double exp = -16.39244950720277;
    const double err = 1e-4;
    EXPECT_LT(abs((res-exp)/exp),err)<<"res = "<<res<<"\t exp = "<<exp<<endl;
}


TEST_F(TestAgainstMathematica,gluon_gluon_reg_L0)
{
    const double z=0.345;
    HiggsMtExpansion_gg_nnlo_reg mtggreg(rho);
    
    double res=mtggreg.z_times_reg_L0(z)/z;
    
    
    const double exp = -274.3984043605598;
    const double err = 1e-4;
    EXPECT_LT(abs((res-exp)/exp),err)<<"res = "<<res<<"\t exp = "<<exp<<endl;
}

TEST_F(TestAgainstMathematica,gluon_gluon_reg_L1)
{
    const double z=0.345;
    HiggsMtExpansion_gg_nnlo_reg mtggreg(rho);
    
    double res=mtggreg.z_times_reg_L1(z)/z;
    
    
    const double exp = -305.6251450939468;
    const double err = 1e-4;
    EXPECT_LT(abs((res-exp)/exp),err)<<"res = "<<res<<"\t exp = "<<exp<<endl;
}




TEST_F(TestAgainstMathematica,gluon_gluon_reg_L2)
{
    const double z=0.345;
    HiggsMtExpansion_gg_nnlo_reg mtggreg(rho);
    
    double res=mtggreg.z_times_reg_L2(z)/z;
    
    
    const double exp = -66.86689023762945;
    const double err = 1e-4;
    EXPECT_LT(abs((res-exp)/exp),err)<<"res = "<<res<<"\t exp = "<<exp<<endl;
}



TEST_F(TestAgainstMathematica,gluon_gluon_delta)
{
   
    
    double res=MTEXP::gg_nnlo_delta(rho,lmh,lmtop);
    
    
    const double exp = 18.51460470739025;
    const double err = 1e-14;
    EXPECT_LT(abs((res-exp)/exp),err)<<"res = "<<res<<"\t exp = "<<exp<<endl;
}


TEST_F(TestAgainstMathematica,gluon_gluon_D0)
{
    
    
    double res=MTEXP::gg_nnlo_D0(rho,lmh);
    
    
    const double exp = 115.113830387372;
    const double err = 1e-14;
    EXPECT_LT(abs((res-exp)/exp),err)<<"res = "<<res<<"\t exp = "<<exp<<endl;
}

TEST_F(TestAgainstMathematica,gluon_gluon_D1)
{
    
    
    double res=MTEXP::gg_nnlo_D1(rho,lmh);
    
    
    const double exp = 269.0099643390618;
    const double err = 1e-14;
    EXPECT_LT(abs((res-exp)/exp),err)<<"res = "<<res<<"\t exp = "<<exp<<endl;
}

TEST_F(TestAgainstMathematica,gluon_gluon_D2)
{
    
    
    double res=MTEXP::gg_nnlo_D2(rho,lmh);
    
    
    const double exp = 324.6385890857657;
    const double err = 1e-14;
    EXPECT_LT(abs((res-exp)/exp),err)<<"res = "<<res<<"\t exp = "<<exp<<endl;
}

TEST_F(TestAgainstMathematica,gluon_gluon_D3)
{
    
    
    double res=MTEXP::gg_nnlo_D3(rho,lmh);
    
    
    const double exp = 72.;
    const double err = 1e-14;
    EXPECT_LT(abs((res-exp)/exp),err)<<"res = "<<res<<"\t exp = "<<exp<<endl;
}






