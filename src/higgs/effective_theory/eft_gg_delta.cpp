
#include "constants.h"
#include<vector>
#include "as_series.h"
#include "higgs_eft.h"
#include <complex>
using namespace std;

// gluon - gluon delta
namespace HEFT{

AsSeries n_delta_at_mh()
{
    

    
    
    return AsSeries(0,
                    n_LO_delta(),
                    n_NLO_delta(),
                    n_NNLO_delta(),
                    n_N3LO_delta()
                    );
}

AsSeries n_delta_log_muf(const double& _log_muf_mh_sq)
{
    return AsSeries(2,
                    n_NNLO_delta_L()*_log_muf_mh_sq
                    + n_NNLO_delta_L2()*pow(_log_muf_mh_sq,2.),
                    n_N3LO_delta_L()*_log_muf_mh_sq
                    + n_N3LO_delta_L2()*pow(_log_muf_mh_sq,2.)
                    + n_N3LO_delta_L3()*pow(_log_muf_mh_sq,3.)
                    );
}

double n_LO_delta(){return 1.0;}

double n_NLO_delta(){return 6.0*consts::z2;}

double n_NNLO_delta(){return
    1071.0 / 8.0 * consts::z4
    - 54.0 * pow(consts::z2,2.0)
    + 67.0/2.0 * consts::z2
    - 165./4. * consts::z3
    + 837./16.
    + ( - 5./3. * consts::z2
       + 5./6.*consts::z3
       - 247./36.)
    * consts::nf;
}
double n_NNLO_delta_L(){return
    33./2. * consts::z2
    -171./2.*consts::z3
    +27./2.
    + (-consts::z2 - 11./6.)
    * consts::nf;
}
double n_NNLO_delta_L2(){return
    - 18. * consts::z2;
}

double n_N3LO_delta(){return
    (-35.*pow(consts::nf,2.)*(-103753.0
                              + 25776.*consts::z2
                              + 5616.*consts::z3
                              + 25200.*consts::z4)
     + 105.*consts::nf*(-1128767.
                        + 429096.*consts::z3
                        - 432.*consts::z2*(188. + 981.*consts::z3)
                        + 287100.*consts::z4
                        + 568872.*consts::z5
                        )
     - 27.*(-420.*consts::z2*(16151. + 52866.*consts::z3)
            + 9611910.*consts::z4
            - 105.*(215131.
                    - 265356.*consts::z3
                    + 356832.*pow(consts::z3,2.)
                    - 272844.*consts::z5)
            + 22714020.*consts::z6)
     )/544320.;
}




    //: numerically equal to -316.3371547476737 (with our convention for the log: L=log(muf^2/mh^2))
    //: value taken from the soft analytic result by Bernhard (November 2017)
double n_N3LO_delta_L(){
    return  (1201. / 576.
             + 5. / 54. * consts::pi_square
             - 5. / 18. * consts::z3) * pow(consts::nf,2.)
    + (-3. / 8. * consts::pi_square
       + 395. / 6. * consts::z3
       - 29807. / 576.
       - 5. / 96. * pow(consts::pi_square,2.)) * consts::nf
    + 18217. / 64.
    + 1089. / 4. * consts::z3 * consts::pi_square
    - 5049. / 2. * consts::z5
    - 11.0 * consts::beta_zero
    + 6. * consts::beta_two
    - 46. / 3. * consts::pi_square
    + 5. * consts::pi_square * consts::beta_one
    - 9453. / 8. * consts::z3
    + 55. / 64. * pow(consts::pi_square,2.)  ;
}
    

    


double n_N3LO_delta_L2(){
    return (-2. / 9.
            + consts::pi_square / 36.)
    * pow(consts::nf,2.)
    + (81. / 4. * consts::z3
       + 17. / 3.
       - 2. / 3. * consts::beta_one
       + 3. / 4.* consts::pi_square)
    * consts::nf
    - 2673. / 8. * consts::z3
    - 27. / 10. * pow(consts::pi_square,2.)
    - 33.
    - 415. / 16. * consts::pi_square
    + 11. * consts::beta_one ;
}

double n_N3LO_delta_L3(){
    return -72. * consts::z3
    - 33. / 4. * consts::pi_square
    + consts::nf * consts::pi_square / 2. ;
}

}
