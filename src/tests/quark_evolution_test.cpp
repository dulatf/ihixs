/** testing UserInterface.*:
 *
 * Achilleas Lazopoulos, lazopoli@phys.ethz.ch
 */

#include <iostream>
#include <cmath>



using namespace std;

#include "gtest/gtest.h"
#include "model.h"
#include "math.h"
class ModelTest: public ::testing::Test {
protected:
    virtual void SetUp()
    {
        as_mz = 0.1170699000000000;
        porder=2; // maximum order is 2
        mur = 125.;
        mh = 125.;
        UserInterface UI;
        UI.SetOption("mb_msbar", "4.163");
        UI.SetOption("mb_msbar_ref_scale", "4.163");
        UI.SetOption("mc_msbar", "0.986");
        UI.SetOption("mc_msbar_ref_scale", "3.0");
        model.ReadParameters(UI);
        model.Configure(as_mz,mur/mh,porder,mh);
    }
    
    void evolve_to(const double& mu){
        mur=mu;
        model.Configure(as_mz,mur/mh,porder,mh);
    }
    double mur;
    double mh;
    double as_mz;
    int porder;
    CModel model;
    
};




TEST_F(ModelTest,mb_at_125)
{
    const double res = model.bottom.m();
    const double exp = 2.77764;
    const double err = 6. * 1e-2;
    EXPECT_LT(abs((res-exp)/exp),err)
    << "our mb(125) = " << res << " vs the nominal value = " << exp;
}

TEST_F(ModelTest,mc_at_125)
{
    const double res = model.charm.m();
    const double exp = 0.620885;
    const double err = 1e-6;
    EXPECT_LT(abs((res-exp)/exp),err)
    << "our mc(125) = " << res << " vs the nominal value = " << exp;
    ;
}


class ModelTestDifferentAsAtMz: public ::testing::Test {
protected:
    virtual void SetUp()
    {
        as_mz = 0.1189;// this is the value they adopt in http://arxiv.org/pdf/0907.2110.pdf to arrive at the mb(161.8) value they quote (which is what we test below)
        porder=2; // maximum order is 2
        mur = 125.;
        mh = 125.;
        UserInterface UI;
        UI.SetOption("mb_msbar", "4.163");
        UI.SetOption("mb_msbar_ref_scale", "4.163");
        UI.SetOption("mc_msbar", "0.986");
        UI.SetOption("mc_msbar_ref_scale", "3.0");
        model.ReadParameters(UI);
        model.Configure(as_mz,mur/mh,porder,mh);
    }
    
    void evolve_to(const double& mu){
        mur=mu;
        model.Configure(as_mz,mur/mh,porder,mh);
    }
    double mur;
    double mh;
    double as_mz;
    int porder;
    CModel model;
    
};

TEST_F(ModelTestDifferentAsAtMz,DISABLED_mb_at_161p8)
{
    evolve_to(161.8);
    const double res = model.bottom.m();
    const double exp = 2.79303;
    const double err = 1e-3;
    EXPECT_LT(abs((res-exp)/exp),err)<<"Failing test indicates that we don't do 6 flavor matching like in http://arxiv.org/pdf/0907.2110.pdf  at the top threshold."<<endl;
}

TEST_F(ModelTestDifferentAsAtMz,mb_at_mz)
{
    evolve_to(91.1876);// mz=91.1876
    const double res = model.bottom.m();
    const double exp = 2.835;
    const double err = 1e-3;
    EXPECT_LT(abs((res-exp)/exp),err)
    << "our mb(91.1876) = " << res << " vs the nominal value = " << exp;
;
}



int main(int argc, char**argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return  RUN_ALL_TESTS();
    return 0;
}




















