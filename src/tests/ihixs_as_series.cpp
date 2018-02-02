/** testing UserInterface.*:
 *
 * Achilleas Lazopoulos, lazopoli@phys.ethz.ch
 */

#include <iostream>
#include <cmath>



using namespace std;

#include "gtest/gtest.h"
#include "as_series.h"
#include "math.h"

class AsSeriesDefaultConstructor: public ::testing::Test {
protected:
    virtual void SetUp()
    {
    }
    
    AsSeries zero;
};


TEST_F(AsSeriesDefaultConstructor,StartingExponent){EXPECT_EQ(zero.starting_exponent(),0);}
TEST_F(AsSeriesDefaultConstructor,EndingExponent){EXPECT_EQ(zero.ending_exponent(),0);}
TEST_F(AsSeriesDefaultConstructor,TermOfOrder0){EXPECT_EQ(zero.term_of_order(0).val(),0.0);}
TEST_F(AsSeriesDefaultConstructor,TermOfOrder1){EXPECT_EQ(zero.term_of_order(1).val(),0.0);}

TEST_F(AsSeriesDefaultConstructor,MulC)
{
    ResultPair c(2.0,0.);
    AsSeries newzero = c*zero;
    EXPECT_EQ(newzero.term_of_order(0).val(),0.0);
}

class SingleMonomial: public ::testing::Test {
protected:
    virtual void SetUp(){pol = AsSeries(1,1.23);}
    AsSeries pol;
};

TEST_F(SingleMonomial,StartingExponent){EXPECT_EQ(pol.starting_exponent(),1);}
TEST_F(SingleMonomial,EndingExponent){EXPECT_EQ(pol.ending_exponent(),1);}
TEST_F(SingleMonomial,TermOfOrder0){EXPECT_EQ(pol.term_of_order(0).val(),0.0);}
TEST_F(SingleMonomial,TermOfOrder1){EXPECT_EQ(pol.term_of_order(1).val(),1.23);}
TEST_F(SingleMonomial,TermOfOrder2){EXPECT_EQ(pol.term_of_order(2).val(),0.);}
TEST_F(SingleMonomial,MulC)
{
    ResultPair c(2.0,0.);
    AsSeries newpol = pol*c;
    EXPECT_EQ(newpol.term_of_order(1).val(),2.46)<<pol<<newpol;
}


TEST_F(SingleMonomial,MulCError)
{
    ResultPair c(2.0,0.);
    AsSeries newpol = c*pol;
    EXPECT_EQ(newpol.term_of_order(1).err(),0.0)<<pol<<newpol;
}

class TwoMonomials: public ::testing::Test {
protected:
    virtual void SetUp(){pol = AsSeries(2,1.23,2.56);}
    AsSeries pol;
};
TEST_F(TwoMonomials,StartingExponent){EXPECT_EQ(pol.starting_exponent(),2);}
TEST_F(TwoMonomials,EndingExponent){EXPECT_EQ(pol.ending_exponent(),3);}
TEST_F(TwoMonomials,TermOfOrder0){EXPECT_EQ(pol.term_of_order(0).val(),0.0);}
TEST_F(TwoMonomials,TermOfOrder1){EXPECT_EQ(pol.term_of_order(1).val(),0.0);}
TEST_F(TwoMonomials,TermOfOrder2){EXPECT_EQ(pol.term_of_order(2).val(),1.23);}
TEST_F(TwoMonomials,TermOfOrder3){EXPECT_EQ(pol.term_of_order(3).val(),2.56);}
TEST_F(TwoMonomials,TermOfOrder4){EXPECT_EQ(pol.term_of_order(4).val(),0.0);}

class ThreeMonomials: public ::testing::Test {
protected:
    virtual void SetUp(){pol = AsSeries(1,1.23,2.56,4.53);}
    AsSeries pol;
};
TEST_F(ThreeMonomials,StartingExponent){EXPECT_EQ(pol.starting_exponent(),1);}
TEST_F(ThreeMonomials,EndingExponent){EXPECT_EQ(pol.ending_exponent(),3);}
TEST_F(ThreeMonomials,TermOfOrder0){EXPECT_EQ(pol.term_of_order(0).val(),0.0);}
TEST_F(ThreeMonomials,TermOfOrder1){EXPECT_EQ(pol.term_of_order(1).val(),1.23);}
TEST_F(ThreeMonomials,TermOfOrder2){EXPECT_EQ(pol.term_of_order(2).val(),2.56);}
TEST_F(ThreeMonomials,TermOfOrder3){EXPECT_EQ(pol.term_of_order(3).val(),4.53);}
TEST_F(ThreeMonomials,TermOfOrder4){EXPECT_EQ(pol.term_of_order(4).val(),0.0);}
TEST_F(ThreeMonomials,MulC)
{
    ResultPair c(2.0,0.);
    AsSeries newpol = pol*c;
    EXPECT_EQ(newpol.term_of_order(1).val(),2.46)<<pol<<newpol;
}
TEST_F(ThreeMonomials,MulC2)
{
    ResultPair c(2.0,0.);
    AsSeries newpol = pol*c;
    EXPECT_EQ(newpol.term_of_order(2).val(),2.56*2)<<pol<<newpol;
}
TEST_F(ThreeMonomials,MulC3)
{
    ResultPair c(2.0,0.);
    AsSeries newpol = pol*c;
    EXPECT_EQ(newpol.term_of_order(3).val(),4.53*2)<<pol<<newpol;
}

TEST_F(ThreeMonomials,MulC4)
{
    ResultPair c(2.0,0.);
    AsSeries newpol = pol*c;
    EXPECT_EQ(newpol.term_of_order(4).val(),0.0)<<pol<<newpol;
}

TEST_F(ThreeMonomials,MulWithDouble)
{
    double c=3.45;
    AsSeries newpol = pol*c;
    EXPECT_EQ(newpol.term_of_order(3).val(),4.53*3.45)<<pol<<newpol;
}

class FourMonomials: public ::testing::Test {
protected:
    virtual void SetUp(){pol = AsSeries(1,1.23,2.56,4.53,75.43);}
    AsSeries pol;
};
TEST_F(FourMonomials,StartingExponent){EXPECT_EQ(pol.starting_exponent(),1);}
TEST_F(FourMonomials,EndingExponent){EXPECT_EQ(pol.ending_exponent(),4);}
TEST_F(FourMonomials,TermOfOrder0){EXPECT_EQ(pol.term_of_order(0).val(),0.0);}
TEST_F(FourMonomials,TermOfOrder1){EXPECT_EQ(pol.term_of_order(1).val(),1.23);}
TEST_F(FourMonomials,TermOfOrder2){EXPECT_EQ(pol.term_of_order(2).val(),2.56);}
TEST_F(FourMonomials,TermOfOrder3){EXPECT_EQ(pol.term_of_order(3).val(),4.53);}
TEST_F(FourMonomials,TermOfOrder4){EXPECT_EQ(pol.term_of_order(4).val(),75.43);}
TEST_F(FourMonomials,TermOfOrder5){EXPECT_EQ(pol.term_of_order(5).val(),0.0);}
TEST_F(FourMonomials,Truncate)
{
    pol.Truncate(2);
    EXPECT_EQ(pol.term_of_order(3).val(),0.0);
}
TEST_F(FourMonomials,Truncate2)
{
    pol.Truncate(2);
    EXPECT_EQ(pol.term_of_order(2).val(),2.56)<<pol;
}
TEST_F(FourMonomials,MultiplyByAs)
{
    pol.MultiplyAs(3.0);
    EXPECT_EQ(pol.term_of_order(2).val(),2.56*pow(3.0,2))<<pol;
}

class AddingSeries: public ::testing::Test {
protected:
    virtual void SetUp(){a = AsSeries(1,1.,2.,3.,4.); b = AsSeries(2,5.,6.);res = a+b;}
    AsSeries a;
    AsSeries b;
    AsSeries res;
};

TEST_F(AddingSeries,Term0){EXPECT_EQ(res.term_of_order(0).val(),0.);}
TEST_F(AddingSeries,Term1){EXPECT_EQ(res.term_of_order(1).val(),1.);}
TEST_F(AddingSeries,Term2){EXPECT_EQ(res.term_of_order(2).val(),7.);}
TEST_F(AddingSeries,Term3){EXPECT_EQ(res.term_of_order(3).val(),9.);}
TEST_F(AddingSeries,Term4){EXPECT_EQ(res.term_of_order(4).val(),4.);}
TEST_F(AddingSeries,Term5){EXPECT_EQ(res.term_of_order(5).val(),0.);}

class MultiplySeries: public ::testing::Test {
protected:
    virtual void SetUp(){a = AsSeries(1,1.,2.,3.,4.); b = AsSeries(2,5.,6.);res = a*b;}
    AsSeries a;
    AsSeries b;
    AsSeries res;
};
TEST_F(MultiplySeries,Term0){EXPECT_EQ(res.term_of_order(0).val(),0.);}
TEST_F(MultiplySeries,Term1){EXPECT_EQ(res.term_of_order(1).val(),0.);}
TEST_F(MultiplySeries,Term2){EXPECT_EQ(res.term_of_order(2).val(),0.);}
TEST_F(MultiplySeries,Term3){EXPECT_EQ(res.term_of_order(3).val(),5.);}
TEST_F(MultiplySeries,Term4){EXPECT_EQ(res.term_of_order(4).val(),16.);}
TEST_F(MultiplySeries,Term5){EXPECT_EQ(res.term_of_order(5).val(),27.);}
TEST_F(MultiplySeries,Term6){EXPECT_EQ(res.term_of_order(6).val(),38.);}
TEST_F(MultiplySeries,Term7){EXPECT_EQ(res.term_of_order(7).val(),24.);}
TEST_F(MultiplySeries,Term8){EXPECT_EQ(res.term_of_order(8).val(),0.);}

class SeriesToThe3: public ::testing::Test {
protected:
    virtual void SetUp(){a = AsSeries(1,1.,2.,3.,4.);res = a.ToThe(3);}
    AsSeries a;
    AsSeries res;
};
TEST_F(SeriesToThe3,Term0){EXPECT_EQ(res.term_of_order(0).val(),0.);}
TEST_F(SeriesToThe3,Term1){EXPECT_EQ(res.term_of_order(1).val(),0.);}
TEST_F(SeriesToThe3,Term2){EXPECT_EQ(res.term_of_order(2).val(),0.);}
TEST_F(SeriesToThe3,Term3){EXPECT_EQ(res.term_of_order(3).val(),1.);}
TEST_F(SeriesToThe3,Term4){EXPECT_EQ(res.term_of_order(4).val(),6.);}
TEST_F(SeriesToThe3,Term5){EXPECT_EQ(res.term_of_order(5).val(),21.);}
TEST_F(SeriesToThe3,Term6){EXPECT_EQ(res.term_of_order(6).val(),56.);}
TEST_F(SeriesToThe3,Term7){EXPECT_EQ(res.term_of_order(7).val(),111.);}
TEST_F(SeriesToThe3,Term8){EXPECT_EQ(res.term_of_order(8).val(),174.);}
TEST_F(SeriesToThe3,Term9){EXPECT_EQ(res.term_of_order(9).val(),219.);}
TEST_F(SeriesToThe3,Term10){EXPECT_EQ(res.term_of_order(10).val(),204.);}
TEST_F(SeriesToThe3,Term11){EXPECT_EQ(res.term_of_order(11).val(),144.);}
TEST_F(SeriesToThe3,Term12){EXPECT_EQ(res.term_of_order(12).val(),64.);}
TEST_F(SeriesToThe3,Term13){EXPECT_EQ(res.term_of_order(13).val(),0.);}


int main(int argc, char**argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return  RUN_ALL_TESTS();
    return 0;
}




















