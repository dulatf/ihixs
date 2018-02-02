/** testing UserInterface.*:
 *
 * Achilleas Lazopoulos, lazopoli@phys.ethz.ch
 */

#include <iostream>
#include <cmath>




using namespace std;

#include "gtest/gtest.h"
#include "nlo_exact_matrix_elements.h"
#include "model.h"
#include <complex>
using namespace std;
TEST(ggf_nlo_exact_virtual,mt_infinity)
{
    complex<double> res,expected,x,y,mq;
    //: if you don't want to see the limiting behavior
    //: comment out from here
    cout<<"\n======================================================"
    <<"\n\t The limit mt->infinity should be 11/2"
    <<"\n\t note that 11/2 is approached well until ~3TeV."
    <<" and then starts deviating, due to cancelations between HPLs"
    <<endl;
    double mh=125.0;
    double mqs[20]={100.0,170.0,300.0,500.0,1000.0,
                    2000.0,2500.0,2800.0,3000.0,3100.0,
                    3200.0,3300.0,3500.0,4000.0,5000.0,
                    10000.0,20000.0,50000.0,75000.0,100000.0};
    for (int i=0;i<20;i++)
        {
        mq= complex<double>(mqs[i],1e-15);

        y =  pow(mq,2.0) / pow(mh,2.0);
        x= ( sqrt(1.0-4.0*y) - 1.0 ) / ( sqrt(1.0-4.0*y) + 1.0 );
        res= 2.0*h_exact::ggf_exact_virtual_ep0(x, 4.0/3.0); //on-shell scheme
        // the limit is scheme-independent
        expected=complex<double>(11.0/2.0,0.0);
    
        cout<<setprecision(16)<<endl
            <<"\tx="<<x<<"\t|x|="<<abs(x)
            <<setprecision(6)<<"\tres="<<real(res)<<"\t|RE(res-expected)|="
            <<setprecision(1)<<scientific<<abs(real(res-expected))
            <<fixed<<"\tmq="<<mq;
        
         }
    cout<<"\n======================================================";
    cout<<endl;
    //: until here
    mq=complex<double>(3200.0,0.0);
    y =  pow(mq,2.0) / pow(mh,2.0);
    x= ( sqrt(1.0-4.0*y) - 1.0 ) / ( sqrt(1.0-4.0*y) + 1.0 );
    res= 2.0*h_exact::ggf_exact_virtual_ep0(x,4.0/3.0); // on-shell scheme
    // the limit is scheme-independent
    expected=complex<double>(11.0/2.0,0.0);
    EXPECT_LT(abs(real(res-expected)),1e-3);
}

TEST(ggf_nlo_exact_virtual,mt_zero)
{
    complex<double> res,expected,x,y,mq;
    //: if you don't want to see the limiting behavior
    //: comment out from here
    cout<<"\n======================================================"
    <<"\n\t The limit mt->0 should be zero"
    <<endl;
    double mh=125.0;
    double mqs[10]={100.0,10.0,5.0,1.0,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6};
    for (int i=0;i<10;i++)
        {
        mq= complex<double>(mqs[i],1e-15);
        
        y =  pow(mq,2.0) / pow(mh,2.0);
        x= ( sqrt(1.0-4.0*y) - 1.0 ) / ( sqrt(1.0-4.0*y) + 1.0 );
        res= 2.0* h_exact::ggf_exact_virtual_ep0(x,4.0/3.0);// on-shell scheme
        // the limit is scheme-independent
        expected=complex<double>(0.0,0.0);
        
        cout<<setprecision(16)<<endl
        <<"\tx="<<x<<"\t|x|="<<abs(x)
        <<setprecision(6)<<scientific<<"\tres="<<real(res)<<"\t|RE(res-expected)|="
        <<setprecision(1)<<scientific<<abs(real(res-expected))
        <<"\tmq="<<mq;
        
        }
    cout<<"\n======================================================";
    cout<<endl;
    //: until here
    mq=complex<double>(1e-4,0.0);
    y =  pow(mq,2.0) / pow(mh,2.0);
    x= ( sqrt(1.0-4.0*y) - 1.0 ) / ( sqrt(1.0-4.0*y) + 1.0 );
    res= 2.0*h_exact::ggf_exact_virtual_ep0(x, 4.0/3.0); // on-shell scheme
    // the limit is scheme-independent
    expected=complex<double>(0.0,0.0);
    EXPECT_LT(abs(real(res-expected)),1e-6);
}

TEST(ggf_nlo_exact_virtual,threshold)
{
    complex<double> res,expected,x,y,mq;
    // threshold limits
    double pisq=consts::pi_square;
    double z3=consts::z3;
    
    double limF2la = -9.0 + 7.0/4.0*pisq - 2.0*pisq*log(2.0) + 7.0*z3;
    double limF2lb = 3.0 - 3.0/4.0*pisq;
    double limG2la = -3.0 + 1.0/2.0*pisq*log(2) - 7.0/4.0*z3;
    expected = complex<double>(-2.0*limF2la -2.0*limF2lb*4.0/3.0
                                - 9.0/2.0 * limG2la,0.0);
    
    
    //: if you don't want to see the limiting behavior
    //: comment out from here
    cout<<"\n======================================================"
    <<"\n\t The limit mt->mh/2 should be "
    <<endl;
    double mh=125.0;
    double mqs[20]={-30.0,-15.0,-1.0,-0.1,-0.01,-0.001,-1e-4,-1e-5,-1e-6,-1e-15,
                    1e-15,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1.0,15.0,30.0};
    for (int i=0;i<20;i++)
        {
        mq= complex<double>(mh/2.0+mqs[i],0.0 );
        
        y =  pow(mq,2.0) / pow(mh,2.0);
        x= ( sqrt(1.0-4.0*y) - 1.0 ) / ( sqrt(1.0-4.0*y) + 1.0 );
        res= 2.0*h_exact::ggf_exact_virtual_ep0(x,4.0/3.0);//on-shell scheme
        
        cout<<setprecision(16)<<endl
        <<"\tx="<<x<<"\t|x|="<<abs(x)
        <<setprecision(6)<<scientific<<"\tres="
        <<real(res)<<"\t|RE(res-expected)|="
        <<setprecision(1)<<scientific<<abs(real(res-expected))
        <<fixed<<setprecision(16)<<"\tmq="<<mq;
        
        }
    cout<<"\n======================================================";
    cout<<"\n---- data for Harlander's plot, fig. 5a in hep-ph/0509189"
    <<"\n\t tau=mh^2/(4*mt^2)"
    <<endl;

    for (int i=0;i<100;i++)
        {

        complex<double> tau =complex<double>( 2.0*double(i)/double(100)+1e-10,0.0);
        y = 1.0/4.0/tau;
        x= ( sqrt(1.0-4.0*y) - 1.0 ) / ( sqrt(1.0-4.0*y) + 1.0 );
        res= 2.0* h_exact::ggf_exact_virtual_ep0(x, 4.0/3.0); // on-shell scheme
        
        cout<<setprecision(16)<<endl
        <<real(1.0/y/4.0)<<" "<<real(res);
        
        }
    cout<<"\n======================================================";
    cout<<endl;
    //: until here
    mq=complex<double>(mh/2.0,0.0);
    y =  pow(mq,2.0) / pow(mh,2.0);
    x= ( sqrt(1.0-4.0*y) - 1.0 ) / ( sqrt(1.0-4.0*y) + 1.0 );
    res= 2.0*h_exact::ggf_exact_virtual_ep0(x, 4.0/3.0); //: on-shell scheme
    cout<<"\n|RE(res-expected)|="<<abs(real(res-expected))<<endl;
    EXPECT_LT(abs(real(res-expected)),1e-6);
}


TEST(ggf_nlo_exact_virtual,threshold_F2lb)
{
    complex<double> res,expected,x,y,mq;
    // threshold limits
    double pisq=consts::pi_square;
    //double z3=consts::z3;
    
    //double limF2la = -9.0 + 7.0/4.0*pisq - 2.0*pisq*log(2.0) + 7.0*z3;
    double limF2lb = 3.0 - 3.0/4.0*pisq;
    //double limG2la = -3.0 + 1.0/2.0*pisq*log(2) - 7.0/4.0*z3;

    double mh=125.0;

    mq=complex<double>(mh/2.0,0.0);
    y =  pow(mq,2.0) / pow(mh,2.0);
    x= ( sqrt(1.0-4.0*y) - 1.0 ) / ( sqrt(1.0-4.0*y) + 1.0 );
    res= h_exact::F2lb(x);
    EXPECT_LT(abs(real(res-limF2lb)),1e-5);
}

class NLOREALTEST: public ::testing::Test {
protected:
    virtual void SetUp()
    {
        Model = new CModel();
        Model->ReadParameters(UI);
        err=1e-14;
        
    }
    
    CModel* Model;
    UserInterface UI;
    double err;
};






TEST_F(NLOREALTEST, ggf_nlo_exact_real_limit_z_to_1)
{
    double lambda = 0.35;
    
    Model->quarks[1]->set_Y(0.0);
    
    Model->Configure(0.117, 0.5, 1,125.0);

    //cout << " hello " << endl;
    
    cout<<"\n---------\n"<<"mh = "<<Model->higgs.m()<<endl;
    cout<<"Contributing quarks: "<<endl;
    for (int i=0;i<Model->quarks.size();i++)
        {
        if (Model->quarks[i]->Y() != 0.0)
            {
            Particle* quark = Model->quarks[i];
            cout<<setprecision(8)
                <<quark->name()<<" "
                <<quark->Y()<<" "
                <<sqrt(quark->cm_sq())
                <<endl;
            }
        }
    cout<<"-----";
    double res,expected_limit;
    
    complex<double> born = h_exact::born_exact_summed_over_quarks(*Model);
    cout<<"\nBorn = "<<born<<endl;
    for (int i=1;i<10;i++)
        {
        double z = 1.0 - 0.1*pow(0.1,double(i));
        
        res = h_exact::sum_of_abs_sq_of_Aqi(z,lambda,
                                            *Model);
        expected_limit = 2.0/pow(z,4.0)*
                                pow(1.0-z+z*z,2.0)*pow(abs(born),2.0);
        cout<<"\n evaluated = "<< res << " expected = "<<expected_limit<<endl;
        }
    EXPECT_LT(abs(res-expected_limit),1e-5);
}



TEST_F(NLOREALTEST, ggf_nlo_exact_real_limit_lambda_to_1)
{
    double z = 0.75;
    Model->quarks[1]->set_Y(0.0);

    Model->Configure(0.117, 0.5, 1,125.0);

    
    cout<<"\n---------\n"<<"mh = "<<Model->higgs.m()<<endl;
    cout<<"Contributing quarks: "<<endl;
    for (int i=0;i<Model->quarks.size();i++)
        {
        if (Model->quarks[i]->Y() != 0.0)
            {
            Particle* quark = Model->quarks[i];
            cout<<setprecision(8)
            <<quark->name()<<" "
            <<quark->Y()<<" "
            <<sqrt(quark->cm_sq())
            <<endl;
            }
        }
    cout<<"-----";
    double res,expected_limit;
    
    complex<double> born = h_exact::born_exact_summed_over_quarks(*Model);
    cout<<"\nBorn = "<<born<<endl;
    for (int i=1;i<8;i++)
        {
        double lambda = 1.0 - 0.1*pow(0.1,double(i));
        
        res = h_exact::sum_of_abs_sq_of_Aqi(z,lambda,
                                            *Model);
        expected_limit = 2.0/pow(z,4.0)*
        pow(1.0-z+z*z,2.0)*pow(abs(born),2.0);
        cout<<"\nlambda="<<lambda
            <<": evaluated = "<< res
            << " expected = "<<expected_limit<<endl;
        }
    EXPECT_LT(abs(res-expected_limit),1e-5);
}

TEST_F(NLOREALTEST, ggf_nlo_exact_real_limit_lambda_to_0)
{
    double z = 0.655;
    Model->quarks[1]->set_Y(0.0);
    Model->Configure(0.117, 0.5, 1,125.0);

    
    cout<<"\n---------\n"<<"mh = "<<Model->higgs.m()<<endl;
    cout<<"Contributing quarks: "<<endl;
    for (int i=0;i<Model->quarks.size();i++)
        {
        if (Model->quarks[i]->Y() != 0.0)
            {
            Particle* quark = Model->quarks[i];
            cout<<setprecision(8)
            <<quark->name()<<" "
            <<quark->Y()<<" "
            <<sqrt(quark->cm_sq())
            <<endl;
            }
        }
    cout<<"-----";
    double res,expected_limit;
    
    complex<double> born = h_exact::born_exact_summed_over_quarks(*Model);
    cout<<"\nBorn = "<<born<<endl;
    for (int i=1;i<8;i++)
        {
        double lambda =  0.1*pow(0.1,double(i));
        
        res = h_exact::sum_of_abs_sq_of_Aqi(z,lambda,
                                            *Model);
        expected_limit = 2.0/pow(z,4.0)*
        pow(1.0-z+z*z,2.0)*pow(abs(born),2.0);
        cout<<"\nlambda="<<lambda
        <<": evaluated = "<< res
        << " expected = "<<expected_limit<<endl;
        }
    EXPECT_LT(abs(res-expected_limit),1e-5);
}


//class TestingNloReal: public CoolInt
//{
//public:
//    TestingNloReal():CoolInt(){};
//    double evaluateIntegral(const double xx[],int i_component);
//    void setModel(CModel* Model){_model = Model;}
//private:
//    CModel* _model;
//};
//
//double TestingNloReal::evaluateIntegral(const double xx[],int i_component)
//{
//    double z=0.732;
//    double lambda = xx[0];
//    complex<double> born = h_exact::born_exact_summed_over_quarks(_model);
//
//    return (
//            pow(z,4.0)/2.0*h_exact::sum_of_abs_sq_of_Aqi(z,lambda,_model)
//            - pow(1.0-z+z*z,2.0) * pow(abs(born),2.0)
//            )
//            /lambda/(1.0-lambda)  + 11.0/6.0 * pow(1.0-z,4.0);
//}


//TEST(ggf_nlo_exact_real, DISABLED_large_mt_limit)
//{
//    CModel* Model = new CModel();
//    Model->quarks[0]->set_pole_mass(2000.0);
//    Model->quarks[1]->set_Y(0.0);
//
//    Model->Configure(0.117, 1.0, 1,125.0);
//
//    
//    TestingNloReal my_dude;
//    my_dude.setModel(Model);
//    my_dude.call_vegas();
//    cout<<"\nres="<<my_dude.result()<<endl;
//    EXPECT_LT(abs(my_dude.result()),1e-5);
//}


TEST(Masters,bubf)
{
    complex<double> m(2.97428550666969     ,-1.681076005913941E-017);
    double s1 = -3421.13630838635;
    double s2 = 15625.7868283835;
    //complex<double> x1 = sqrt(1.0-4.0*m*m/s1);
    //complex<double> x2 = sqrt(1.0-4.0*m*m/s2);
    complex<double> expected(  1.47343745569861,-3.13803347841651);
    complex<double> res = h_exact::bubf(m,s1,s2);
    EXPECT_LT(abs(real(res-expected))/abs(expected),1e-15)<<" expected = "<<expected
    <<"\t res = "<<res;
}

TEST(Masters,triaf1)
{
    complex<double> m(2.97428550666969     ,-1.681076005913941E-017);
    double s1 = -1394.30085241407;
    complex<double> expected(  12.8661981549883,-8.095806348395888E-013);
    complex<double> res = h_exact::spec_triaf(sqrt(1.0-4.0*m*m/s1));
    EXPECT_LT(abs(real(res-expected))/abs(expected),1e-15)<<" expected = "<<expected
    <<"\t res = "<<res;
}

TEST(Masters,triaf2)
{
    complex<double> m(2.97428550666969     ,-1.681076005913941E-017);
    double s1 = 15624.5443840588;
    complex<double> expected(  23.0064286786935 , -23.4848417178919);
    complex<double> res = h_exact::spec_triaf(sqrt(1.0-4.0*m*m/s1));
    EXPECT_LT(abs(real(res-expected))/abs(expected),1e-14)<<" expected = "<<expected
    <<"\t res = "<<res;
}

TEST(Aqqgh,bottom_quark)
{
    double z = -7.23453383973682;
    //double mh = 125.003249603365;
    complex<double> tau ( 2.264554070476561E-003,-2.559866901436423E-020);
    complex<double> x = (sqrt(1.0-tau)-1.0)/(sqrt(1.0-tau)+1.0);
    complex<double> m(2.97428550666969     ,-1.681076005913941E-017);

    complex<double> expected(  -4.71636081027699     ,  19.9224142639752);
    complex<double> res = h_exact::Aqqgh_cpp(z,x);
    EXPECT_LT(abs(real(res-expected))/abs(expected),1e-13)<<" expected = "<<expected
    <<"\t res = "<<res;
}


//-706.168901481198      ( 0.996067101904021     , 3.972710522846434E-002)
//-4.88975870452380      (  1.01600912524306     , 3.157339729178606E-002)





/*

TEST_F(NLOREALTEST,sum_Aqqgh_StandardModel)
{
    Model->Configure(0.120180436153365,0.4567,1,125.0);//: a_s@mz, mur/mh, porder
    double z = -706.168901481198;

    complex<double> expected_cplx(0.996067101904021,3.97271052284643E-002);
    double expected = pow(abs(expected_cplx),2.0);
    double res = h_exact::sum_of_abs_sq_of_Aqqgh(z,*Model);
    EXPECT_LT(abs((res-expected))/abs(expected),1e-5)<<" expected = "<<expected
    <<"\t res = "<<res;
}

TEST_F(NLOREALTEST,sum_AqqghStandardModel2)
{
    Model->Configure(0.120180436153365,0.4567,1,125.0);//: a_s@mz, mur/mh, porder
    double z = -4.88975870452380;
    
    complex<double> expected_cplx(1.01600912524306,3.157339729178606E-002);
    double expected = pow(abs(expected_cplx),2.0);
    double res = h_exact::sum_of_abs_sq_of_Aqqgh(z,*Model);
    EXPECT_LT(abs((res-expected))/abs(expected),1e-4)<<" expected = "<<expected
    <<"\t res = "<<res;
}


TEST_F(NLOREALTEST,ggRegNLOExact_LisZero)
{
    // the gg_reg_exact_nlo depends on z, lambda and L
    // but also on the model. In this test we have the default model
    // assuming that this has msbar masses for top and bottom
    // otherwise the test will fail
    const double z=0.23486328125000000;
    const double lambda = 0.11474609375000000;
    const double expected = 20.284340752826505;
    const double L=0;
    const double err=1e-5;
    const double as_mz = 0.11707;
    const double m_h = 125.0;
    const double mur_over_mh = 1.;
    Model->Configure(
                     as_mz,
                     mur_over_mh,
                     2,
                     m_h
                     );
    const double res = h_exact::gg_reg_exact_nlo(z,lambda,L,*Model);
    EXPECT_LT(abs(res-expected)/abs(expected),err)
    <<" res="<<res
    <<" expected="<<expected
    <<"\t % diff = "<<abs(res-expected)/abs(expected)
    <<endl;
}
*/


int main(int argc, char**argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return  RUN_ALL_TESTS();
    return 0;
}
