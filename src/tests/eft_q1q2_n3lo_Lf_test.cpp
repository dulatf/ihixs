/** testing UserInterface.*:
 *
 * Achilleas Lazopoulos, lazopoli@phys.ethz.ch
 */

#include <iostream>
#include <cmath>



using namespace std;

#include "gtest/gtest.h"
#include "higgs_eft.h"
#include "math.h"

// testing a new implementation of LEq1q2N3LOregFalko_L3

TEST(n3Lq1q2simplify,L3)
{
    double maximum_relative_error = 0.0;
    int number_of_points = 100;
    for ( int i=0; i < number_of_points - 1; i++ ){
        double zz=(i+1)/ double(number_of_points);
        
        double exact = HEFT::LEq1q2N3LOregFalko_L3(zz);
        double new_implementation = HEFT::LEq1q2N3LOregFalko_L3(zz);
        
        double relative_error = (exact-new_implementation)/exact ;
        maximum_relative_error = max(relative_error,maximum_relative_error);
        cout << "zz=" << zz << "\t" << new_implementation
        << "\t" << exact
        << "\t delta(\%) = " << relative_error * 100.
        << endl;
    }
    
    
    const double res = maximum_relative_error;
    const double exp = 1e-12;
    cout << "\nprecision achieved : " << maximum_relative_error << endl;
    EXPECT_LT(res,exp) << " res = " << res << " expected = " << exp << endl;
}

TEST(n3Lq1q2simplify,L2)
{
    double maximum_relative_error = 0.0;
    int number_of_points = 100;
    for ( int i=0; i < number_of_points - 1; i++ ){
        double zz=(i+1)/ double(number_of_points);
        
        double exact = HEFT::LEq1q2N3LOregFalko_L2(zz);
        double new_implementation = HEFT::LEq1q2N3LOregFalko_L2(zz);
        
        double relative_error = (exact-new_implementation)/exact ;
        maximum_relative_error = max(relative_error,maximum_relative_error);
        cout << "zz=" << zz << "\t" << new_implementation
        << "\t" << exact
        << "\t delta(\%) = " << relative_error * 100.
        << endl;
    }
    
    
    const double res = maximum_relative_error;
    const double exp = 1e-12;
    cout << "\nprecision achieved : " << maximum_relative_error << endl;
    EXPECT_LT(res,exp) << " res = " << res << " expected = " << exp << endl;
}

TEST(n3Lq1q2simplify,L1)
{
    double maximum_relative_error = 0.0;
    int number_of_points = 100;
    for ( int i=0; i < number_of_points - 1; i++ ){
        double zz=(i+1)/ double(number_of_points);
        
        double exact = HEFT::LEq1q2N3LOregFalko_L1(zz);
        double new_implementation = HEFT::LEq1q2N3LOregFalko_L1(zz);
        
        double relative_error = (exact-new_implementation)/exact ;
        maximum_relative_error = max(relative_error,maximum_relative_error);
        cout << "zz=" << zz << "\t" << new_implementation
        << "\t" << exact
        << "\t delta(\%) = " << relative_error * 100.
        << endl;
    }
    
    
    const double res = maximum_relative_error;
    const double exp = 1e-12;
    cout << "\nprecision achieved : " << maximum_relative_error << endl;
    EXPECT_LT(res,exp) << " res = " << res << " expected = " << exp << endl;
}


int main(int argc, char**argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return  RUN_ALL_TESTS();
    return 0;
}




















