
#include "constants.h"
#include<vector>
#include "as_series.h"
#include "higgs_eft.h"
#include <complex>
using namespace std;

// gluon - gluon plus
namespace HEFT{

    
    
    AsSeries n_D0_at_mh()
    {
        return AsSeries(1,
                        0.,
                        n_NNLO_D0(),
                        n_N3LO_D0()
                        );
    }
    
    AsSeries n_D0_log_muf(const double& _log_muf_mh_sq)
    {
        return AsSeries(1,
                        n_NLO_D0_L() * _log_muf_mh_sq,
                        n_NNLO_D0_L()*_log_muf_mh_sq
                        + n_NNLO_D0_L2()*pow(_log_muf_mh_sq,2.),
                        n_N3LO_D0_L()*_log_muf_mh_sq
                        + n_N3LO_D0_L2()*pow(_log_muf_mh_sq,2.)
                        + n_N3LO_D0_L3()*pow(_log_muf_mh_sq,3.)
                        );
    }
    
    AsSeries n_D1_at_mh()
    {
        return AsSeries(1,
                        n_NLO_D1(),
                        n_NNLO_D1(),
                        n_N3LO_D1()
                        );
    }
    
    AsSeries n_D1_log_muf(const double& _log_muf_mh_sq)
    {
        return AsSeries(1,
                        0.0,
                        n_NNLO_D1_L()*_log_muf_mh_sq
                        + n_NNLO_D1_L2()*pow(_log_muf_mh_sq,2.),
                        n_N3LO_D1_L()*_log_muf_mh_sq
                        + n_N3LO_D1_L2()*pow(_log_muf_mh_sq,2.)
                        + n_N3LO_D1_L3()*pow(_log_muf_mh_sq,3.)
                        );
    }
    
    AsSeries n_D2_at_mh()
    {
        return AsSeries(2,
                        n_NNLO_D2(),
                        n_N3LO_D2()
                        );
    }
    
    AsSeries n_D2_log_muf(const double& _log_muf_mh_sq)
    {
        return AsSeries(2,
                        n_NNLO_D2_L()*_log_muf_mh_sq,
                        n_N3LO_D2_L()*_log_muf_mh_sq
                        + n_N3LO_D2_L2()*pow(_log_muf_mh_sq,2.)
                        + n_N3LO_D2_L3()*pow(_log_muf_mh_sq,3.)
                        );
    }
    
    AsSeries n_D3_at_mh()
    {
        return AsSeries(2,
                        n_NNLO_D3(),
                        n_N3LO_D3()
                        );
    }
    
    AsSeries n_D3_log_muf(const double& _log_muf_mh_sq)
    {
        return AsSeries(2,
                        0.,
                        n_N3LO_D3_L()*_log_muf_mh_sq
                        + n_N3LO_D3_L2()*pow(_log_muf_mh_sq,2.)
                        );
    }
    
    AsSeries n_D4_at_mh()
    {
        return AsSeries(3,
                        n_N3LO_D4());
    }
    
    AsSeries n_D4_log_muf(const double& _log_muf_mh_sq)
    {
        return AsSeries(3,
                        n_N3LO_D4_L()*_log_muf_mh_sq);
    }
    
    AsSeries n_D5_at_mh()
    {
        return AsSeries(3,n_N3LO_D5());
    }
    
    
    
    
    
    double n_NLO_D0_L(){return -6.;}
    
    double n_NLO_D1(){return 12.;}
    
    double n_NNLO_D0(){return   - 2.  * consts::z2 * consts::nf
        + 14. / 9. * consts::nf
        + 351. / 2.  * consts::z3
        + 33.  * consts::z2
        - 101. / 3. ;}
    double n_NNLO_D0_L(){return + 5. / 3. * consts::nf
        + 45.  * consts::z2
        - 67./ 2. ;}
    double n_NNLO_D0_L2(){return  1./ 2. * consts::nf - 33. / 4.; }
    
    double n_NNLO_D1(){return   - 10. / 3. * consts::nf
        - 90.  * consts::z2
        + 67. ;
    }
    double n_NNLO_D1_L(){return - 2. * consts::nf + 33.;}
    double n_NNLO_D1_L2(){return + 36.;}
    
    double n_NNLO_D2(){return - 33. + 2. *  consts::nf ;}
    double n_NNLO_D2_L(){return - 108.;}
    double n_NNLO_D3(){return 72.;}
    
    
    double n_N3LO_D0(){return ( -58./243.+(10./9.)*consts::z2
                               +(5./9.)*consts::z3)*pow(consts::nf,2.)
        +(
          41579./1296.-(2245./36.)*consts::z2
          -(4427./36.)*consts::z3
          -(59./10.)*pow(consts::z2,2.)
          )*consts::nf
        +(
          -297029./864.+(8563./12.)*consts::z2
          +(8941./4.)*consts::z3
          +(2277./20.)*pow(consts::z2,2.)
          -(6525./2.)*consts::z2*consts::z3
          +5022.*consts::z5
          );
    }
    double n_N3LO_D0_L(){return  (consts::pi_square / 9.
                                  - 25. / 54.
                                  )* pow(consts::nf,2.)
        + (  5345. / 72.
           - 81. *consts::z3
           - 47. / 6. * consts::pi_square
           ) * consts::nf
        - 30569. / 48.
        + 114. * consts::pi_square
        + 1584. *consts::z3
        + 231. / 20. * pow(consts::pi_square,2.)
        ;
    }
    double n_N3LO_D0_L2(){return - 5. / 18. * pow(consts::nf,2.)
        + (   13. / 6.
           - 7. / 4. * consts::pi_square
           ) *consts::nf
        + 945. *consts::z3
        + 231. / 8. * consts::pi_square
        - 161. / 8.
        - 27. *consts::beta_one ;
    }
    double n_N3LO_D0_L3(){return - pow(consts::nf,2.)/ 18.
        + 11. / 6. *consts::nf
        + (pow(12.*consts::Pi,2) - 121.) / 8. ;
    }
    
    double n_N3LO_D1(){return   (25./27.-(4./3.)*consts::z2)*pow(consts::nf,2.)
        +(-5345./36.+162.*consts::z3
          +94.*consts::z2
          )*consts::nf
        +30569./24.
        -1368.*consts::z2
        -3168.*consts::z3
        -(4158./5.)*pow(consts::z2,2.)
        ;}
    double n_N3LO_D1_L(){return   10. / 9. * pow(consts::nf,2.)
        + (-130. / 3.+13. * consts::pi_square) *consts::nf
        + 60. *consts::beta_one
        -429. /2. * consts::pi_square
        -4860. *consts::z3
        +1257. /2.
        ;}
    double n_N3LO_D1_L2(){return 1971. /4.
        + pow(consts::nf,2.) /3.
        -162. * consts::pi_square
        -31. *consts::nf ;}
    double n_N3LO_D1_L3(){return -6. *consts::nf +99. ;}
    
    double n_N3LO_D2(){return   -1051.
        +1683.*consts::z2
        +4887.*consts::z3
        +(469./6.-102.*consts::z2)*consts::nf
        -(10./9.)*pow(consts::nf,2.)
        ;}
    double n_N3LO_D2_L(){return  -2775. /2.
        +378. * consts::pi_square
        -2. /3. * pow(consts::nf,2.)
        +82. *consts::nf
        ;}
    double n_N3LO_D2_L2(){return + (-891. /2. +27. * consts::nf);}
    double n_N3LO_D2_L3(){return-108. ;}
    
    double n_N3LO_D3(){return 925.
        -1512.*consts::z2
        -(164./3.)*consts::nf
        +(4./9.)*pow(consts::nf,2.)
        ;}
    double n_N3LO_D3_L(){return + (660. -40. * consts::nf) ;}
    double n_N3LO_D3_L2(){return + 432. ;}
    
    double n_N3LO_D4(){return  -330.+20.*consts::nf;}
    double n_N3LO_D4_L(){return -540. ;}
    
    double n_N3LO_D5(){return 216.;}
    

    
}