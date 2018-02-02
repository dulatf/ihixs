#include <strstream>
#include "threshold_resummation.h"
#include "constants.h"
#include "psi.h"
// allocating the null pointer
//PsiFinder *PsiFinder::s_instance = 0;
using namespace std;

namespace RES{




    complex<double> LogN(const complex<double>& n, const string& resummation_type){
        if (resummation_type=="log"){
            //cout<<real(n)<<"+I*"<<imag(n)<<","<<endl;
            return log(n);
        }
        else if (resummation_type=="psi"){
            return psi(n);
        }
        else if (resummation_type=="AP2log"){
            return 2.*log(n)-3.*log(n+1.)+2.*log(n+2.);
        }
        else if (resummation_type=="AP2psi"){
            return 2. * psi(n)
                  -3. * psi(n+1.)
                  +2. * psi(n+2.);
        }
        else{
            cout<<"unknown resummation flavor "<<resummation_type<<endl;
            exit(0);
        }
    }

  // these are the n^0 terms of the xsections in mellin space
  double m_LO_(void){
    return 1.0;
  }
  double m_NLO_(double Lmh2Omuf2, double Lmuf2Omur2){
    return 2*pow(consts::egamma,2)*QCD::Nc - 2*consts::egamma*Lmh2Omuf2*QCD::Nc + (2*QCD::Nc*pow(consts::Pi,2))/3.;
  }
  double m_NNLO_(double Lmh2Omuf2, double Lmuf2Omur2){
    double l2c = (11*consts::egamma*pow(QCD::Nc,2))/12. + 2*pow(consts::egamma,2)*pow(QCD::Nc,2) - (consts::egamma*QCD::Nc*QCD::Nf)/6.;
    double l1c = (-3*pow(QCD::Nc,2))/2. - (67*consts::egamma*pow(QCD::Nc,2))/18. - (11*pow(consts::egamma,2)*pow(QCD::Nc,2))/6. - 4*pow(consts::egamma,3)*pow(QCD::Nc,2) + (11*consts::egamma*Lmuf2Omur2*pow(QCD::Nc,2))/2. - QCD::Nf/(8.*QCD::Nc) + (5*QCD::Nc*QCD::Nf)/8. + (5*consts::egamma*QCD::Nc*QCD::Nf)/9. + (pow(consts::egamma,2)*QCD::Nc*QCD::Nf)/3. - consts::egamma*Lmuf2Omur2*QCD::Nc*QCD::Nf - (11*pow(QCD::Nc,2)*pow(consts::Pi,2))/18. - (7*consts::egamma*pow(QCD::Nc,2)*pow(consts::Pi,2))/6. + (QCD::Nc*QCD::Nf*pow(consts::Pi,2))/9. + (3*pow(QCD::Nc,2)*consts::z3)/2.;
    double l0c = (93*pow(QCD::Nc,2))/16. + (101*consts::egamma*pow(QCD::Nc,2))/27. + (67*pow(consts::egamma,2)*pow(QCD::Nc,2))/18. + (11*pow(consts::egamma,3)*pow(QCD::Nc,2))/9. + 2*pow(consts::egamma,4)*pow(QCD::Nc,2) - (11*pow(consts::egamma,2)*Lmuf2Omur2*pow(QCD::Nc,2))/2. + (67*QCD::Nf)/(96.*QCD::Nc) - (227*QCD::Nc*QCD::Nf)/96. - (14*consts::egamma*QCD::Nc*QCD::Nf)/27. - (5*pow(consts::egamma,2)*QCD::Nc*QCD::Nf)/9. - (2*pow(consts::egamma,3)*QCD::Nc*QCD::Nf)/9. + pow(consts::egamma,2)*Lmuf2Omur2*QCD::Nc*QCD::Nf + (67*pow(QCD::Nc,2)*pow(consts::Pi,2))/54. + (7*pow(consts::egamma,2)*pow(QCD::Nc,2)*pow(consts::Pi,2))/6. - (11*Lmuf2Omur2*pow(QCD::Nc,2)*pow(consts::Pi,2))/6. - (5*QCD::Nc*QCD::Nf*pow(consts::Pi,2))/27. + (Lmuf2Omur2*QCD::Nc*QCD::Nf*pow(consts::Pi,2))/3. + (23*pow(QCD::Nc,2)*pow(consts::Pi,4))/144. - (77*pow(QCD::Nc,2)*consts::z3)/36. - (7*consts::egamma*pow(QCD::Nc,2)*consts::z3)/2. - (QCD::Nf*consts::z3)/(2.*QCD::Nc) - (QCD::Nc*QCD::Nf*consts::z3)/9.;
    return pow(Lmh2Omuf2,2)*l2c + pow(Lmh2Omuf2,1)*l1c + pow(Lmh2Omuf2,0)*l0c;
  }
  double m_NNNLO_(double Lmh2Omuf2, double Lmuf2Omur2){
    double l3c = (-121*consts::egamma*pow(QCD::Nc,3))/216. - (11*pow(consts::egamma,2)*pow(QCD::Nc,3))/6. - (4*pow(consts::egamma,3)*pow(QCD::Nc,3))/3. + (11*consts::egamma*pow(QCD::Nc,2)*QCD::Nf)/54. + (pow(consts::egamma,2)*pow(QCD::Nc,2)*QCD::Nf)/3. - (consts::egamma*QCD::Nc*pow(QCD::Nf,2))/54.;
    double l2c = (11*pow(QCD::Nc,3))/8. + (769*consts::egamma*pow(QCD::Nc,3))/108. + (73*pow(consts::egamma,2)*pow(QCD::Nc,3))/8. + (11*pow(consts::egamma,3)*pow(QCD::Nc,3))/2. + 4*pow(consts::egamma,4)*pow(QCD::Nc,3) - (121*consts::egamma*Lmuf2Omur2*pow(QCD::Nc,3))/36. - (22*pow(consts::egamma,2)*Lmuf2Omur2*pow(QCD::Nc,3))/3. + (11*QCD::Nf)/96. + (5*consts::egamma*QCD::Nf)/16. - (79*pow(QCD::Nc,2)*QCD::Nf)/96. - (1145*consts::egamma*pow(QCD::Nc,2)*QCD::Nf)/432. - (31*pow(consts::egamma,2)*pow(QCD::Nc,2)*QCD::Nf)/18. - pow(consts::egamma,3)*pow(QCD::Nc,2)*QCD::Nf + (11*consts::egamma*Lmuf2Omur2*pow(QCD::Nc,2)*QCD::Nf)/9. + (4*pow(consts::egamma,2)*Lmuf2Omur2*pow(QCD::Nc,2)*QCD::Nf)/3. - pow(QCD::Nf,2)/(48.*QCD::Nc) + (5*QCD::Nc*pow(QCD::Nf,2))/48. + (5*consts::egamma*QCD::Nc*pow(QCD::Nf,2))/54. + (pow(consts::egamma,2)*QCD::Nc*pow(QCD::Nf,2))/18. - (consts::egamma*Lmuf2Omur2*QCD::Nc*pow(QCD::Nf,2))/9. + (121*pow(QCD::Nc,3)*pow(consts::Pi,2))/216. + (121*consts::egamma*pow(QCD::Nc,3)*pow(consts::Pi,2))/72. + pow(consts::egamma,2)*pow(QCD::Nc,3)*pow(consts::Pi,2) - (11*pow(QCD::Nc,2)*QCD::Nf*pow(consts::Pi,2))/54. - (11*consts::egamma*pow(QCD::Nc,2)*QCD::Nf*pow(consts::Pi,2))/36. + (QCD::Nc*pow(QCD::Nf,2)*pow(consts::Pi,2))/54. - (11*pow(QCD::Nc,3)*consts::z3)/8. - 3*consts::egamma*pow(QCD::Nc,3)*consts::z3 + (pow(QCD::Nc,2)*QCD::Nf*consts::z3)/4.;

    double l1c = (-2071*pow(QCD::Nc,3))/144. - (30569*consts::egamma*pow(QCD::Nc,3))/1296. - (337*pow(consts::egamma,2)*pow(QCD::Nc,3))/18. - (925*pow(consts::egamma,3)*pow(QCD::Nc,3))/54. - (55*pow(consts::egamma,4)*pow(QCD::Nc,3))/9. - 4*pow(consts::egamma,5)*pow(QCD::Nc,3) + (11*Lmuf2Omur2*pow(QCD::Nc,3))/2. + (1933*consts::egamma*Lmuf2Omur2*pow(QCD::Nc,3))/108. + (121*pow(consts::egamma,2)*Lmuf2Omur2*pow(QCD::Nc,3))/18. + (44*pow(consts::egamma,3)*Lmuf2Omur2*pow(QCD::Nc,3))/3. - (121*consts::egamma*pow(Lmuf2Omur2,2)*pow(QCD::Nc,3))/12. - (151*QCD::Nf)/96. - (63*consts::egamma*QCD::Nf)/32. - (3*pow(consts::egamma,2)*QCD::Nf)/8. + (11*Lmuf2Omur2*QCD::Nf)/24. + (3*consts::egamma*Lmuf2Omur2*QCD::Nf)/8. - QCD::Nf/(64.*pow(QCD::Nc,2)) + (4973*pow(QCD::Nc,2)*QCD::Nf)/576. + (21947*consts::egamma*pow(QCD::Nc,2)*QCD::Nf)/2592. + (1099*pow(consts::egamma,2)*pow(QCD::Nc,2)*QCD::Nf)/216. + (82*pow(consts::egamma,3)*pow(QCD::Nc,2)*QCD::Nf)/27. + (10*pow(consts::egamma,4)*pow(QCD::Nc,2)*QCD::Nf)/9. - (79*Lmuf2Omur2*pow(QCD::Nc,2)*QCD::Nf)/24. - (1327*consts::egamma*Lmuf2Omur2*pow(QCD::Nc,2)*QCD::Nf)/216. - (22*pow(consts::egamma,2)*Lmuf2Omur2*pow(QCD::Nc,2)*QCD::Nf)/9. - (8*pow(consts::egamma,3)*Lmuf2Omur2*pow(QCD::Nc,2)*QCD::Nf)/3. + (11*consts::egamma*pow(Lmuf2Omur2,2)*pow(QCD::Nc,2)*QCD::Nf)/3. + (13*pow(QCD::Nf,2))/(48.*QCD::Nc) - (Lmuf2Omur2*pow(QCD::Nf,2))/(12.*QCD::Nc) - (263*QCD::Nc*pow(QCD::Nf,2))/288. - (25*consts::egamma*QCD::Nc*pow(QCD::Nf,2))/162. - (5*pow(consts::egamma,2)*QCD::Nc*pow(QCD::Nf,2))/27. - (2*pow(consts::egamma,3)*QCD::Nc*pow(QCD::Nf,2))/27. + (5*Lmuf2Omur2*QCD::Nc*pow(QCD::Nf,2))/12. + (10*consts::egamma*Lmuf2Omur2*QCD::Nc*pow(QCD::Nf,2))/27. + (2*pow(consts::egamma,2)*Lmuf2Omur2*QCD::Nc*pow(QCD::Nf,2))/9. - (consts::egamma*pow(Lmuf2Omur2,2)*QCD::Nc*pow(QCD::Nf,2))/3. - (2419*pow(QCD::Nc,3)*pow(consts::Pi,2))/648. - (469*consts::egamma*pow(QCD::Nc,3)*pow(consts::Pi,2))/108. - (77*pow(consts::egamma,2)*pow(QCD::Nc,3)*pow(consts::Pi,2))/36. - 2*pow(consts::egamma,3)*pow(QCD::Nc,3)*pow(consts::Pi,2) + (121*Lmuf2Omur2*pow(QCD::Nc,3)*pow(consts::Pi,2))/54. + (77*consts::egamma*Lmuf2Omur2*pow(QCD::Nc,3)*pow(consts::Pi,2))/18. - (QCD::Nf*pow(consts::Pi,2))/8. + (433*pow(QCD::Nc,2)*QCD::Nf*pow(consts::Pi,2))/324. + (35*consts::egamma*pow(QCD::Nc,2)*QCD::Nf*pow(consts::Pi,2))/54. + (7*pow(consts::egamma,2)*pow(QCD::Nc,2)*QCD::Nf*pow(consts::Pi,2))/18. - (22*Lmuf2Omur2*pow(QCD::Nc,2)*QCD::Nf*pow(consts::Pi,2))/27. - (7*consts::egamma*Lmuf2Omur2*pow(QCD::Nc,2)*QCD::Nf*pow(consts::Pi,2))/9. - (5*QCD::Nc*pow(QCD::Nf,2)*pow(consts::Pi,2))/81. + (2*Lmuf2Omur2*QCD::Nc*pow(QCD::Nf,2)*pow(consts::Pi,2))/27. - (55*pow(QCD::Nc,3)*pow(consts::Pi,4))/192. - (43*consts::egamma*pow(QCD::Nc,3)*pow(consts::Pi,4))/180. + (5*pow(QCD::Nc,2)*QCD::Nf*pow(consts::Pi,4))/96. + (2053*pow(QCD::Nc,3)*consts::z3)/216. + (88*consts::egamma*pow(QCD::Nc,3)*consts::z3)/9. + 10*pow(consts::egamma,2)*pow(QCD::Nc,3)*consts::z3 - (11*Lmuf2Omur2*pow(QCD::Nc,3)*consts::z3)/2. + (11*QCD::Nf*consts::z3)/12. + (3*consts::egamma*QCD::Nf*consts::z3)/2. - (145*pow(QCD::Nc,2)*QCD::Nf*consts::z3)/108. - (5*consts::egamma*pow(QCD::Nc,2)*QCD::Nf*consts::z3)/18. + Lmuf2Omur2*pow(QCD::Nc,2)*QCD::Nf*consts::z3 - (pow(QCD::Nf,2)*consts::z3)/(6.*QCD::Nc) - (QCD::Nc*pow(QCD::Nf,2)*consts::z3)/27. + (11*pow(QCD::Nc,3)*pow(consts::Pi,2)*consts::z3)/12. - (5*pow(QCD::Nc,3)*consts::z5)/2.;

    double l0c = (215131*pow(QCD::Nc,3))/5184. + (297029*consts::egamma*pow(QCD::Nc,3))/23328. + (30569*pow(consts::egamma,2)*pow(QCD::Nc,3))/1296. + (1051*pow(consts::egamma,3)*pow(QCD::Nc,3))/81. + (925*pow(consts::egamma,4)*pow(QCD::Nc,3))/108. + (22*pow(consts::egamma,5)*pow(QCD::Nc,3))/9. + (4*pow(consts::egamma,6)*pow(QCD::Nc,3))/3. - (341*Lmuf2Omur2*pow(QCD::Nc,3))/16. - (1111*consts::egamma*Lmuf2Omur2*pow(QCD::Nc,3))/81. - (1933*pow(consts::egamma,2)*Lmuf2Omur2*pow(QCD::Nc,3))/108. - (121*pow(consts::egamma,3)*Lmuf2Omur2*pow(QCD::Nc,3))/27. - (22*pow(consts::egamma,4)*Lmuf2Omur2*pow(QCD::Nc,3))/3. + (121*pow(consts::egamma,2)*pow(Lmuf2Omur2,2)*pow(QCD::Nc,3))/12. + (58519*QCD::Nf)/10368. + (1711*consts::egamma*QCD::Nf)/1728. + (63*pow(consts::egamma,2)*QCD::Nf)/32. + (pow(consts::egamma,3)*QCD::Nf)/12. - (737*Lmuf2Omur2*QCD::Nf)/288. - (3*pow(consts::egamma,2)*Lmuf2Omur2*QCD::Nf)/8. + (19*QCD::Nf)/(72.*pow(QCD::Nc,2)) - (28597*pow(QCD::Nc,2)*QCD::Nf)/1152. - (171449*consts::egamma*pow(QCD::Nc,2)*QCD::Nf)/46656. - (21947*pow(consts::egamma,2)*pow(QCD::Nc,2)*QCD::Nf)/2592. - (941*pow(consts::egamma,3)*pow(QCD::Nc,2)*QCD::Nf)/324. - (41*pow(consts::egamma,4)*pow(QCD::Nc,2)*QCD::Nf)/27. - (4*pow(consts::egamma,5)*pow(QCD::Nc,2)*QCD::Nf)/9. + (3613*Lmuf2Omur2*pow(QCD::Nc,2)*QCD::Nf)/288. + (356*consts::egamma*Lmuf2Omur2*pow(QCD::Nc,2)*QCD::Nf)/81. + (1327*pow(consts::egamma,2)*Lmuf2Omur2*pow(QCD::Nc,2)*QCD::Nf)/216. + (44*pow(consts::egamma,3)*Lmuf2Omur2*pow(QCD::Nc,2)*QCD::Nf)/27. + (4*pow(consts::egamma,4)*Lmuf2Omur2*pow(QCD::Nc,2)*QCD::Nf)/3. - (11*pow(consts::egamma,2)*pow(Lmuf2Omur2,2)*pow(QCD::Nc,2)*QCD::Nf)/3. - (4481*pow(QCD::Nf,2))/(5184.*QCD::Nc) + (67*Lmuf2Omur2*pow(QCD::Nf,2))/(144.*QCD::Nc) + (6013*QCD::Nc*pow(QCD::Nf,2))/2592. + (58*consts::egamma*QCD::Nc*pow(QCD::Nf,2))/729. + (25*pow(consts::egamma,2)*QCD::Nc*pow(QCD::Nf,2))/162. + (10*pow(consts::egamma,3)*QCD::Nc*pow(QCD::Nf,2))/81. + (pow(consts::egamma,4)*QCD::Nc*pow(QCD::Nf,2))/27. - (227*Lmuf2Omur2*QCD::Nc*pow(QCD::Nf,2))/144. - (28*consts::egamma*Lmuf2Omur2*QCD::Nc*pow(QCD::Nf,2))/81. - (10*pow(consts::egamma,2)*Lmuf2Omur2*QCD::Nc*pow(QCD::Nf,2))/27. - (4*pow(consts::egamma,3)*Lmuf2Omur2*QCD::Nc*pow(QCD::Nf,2))/27. + (pow(consts::egamma,2)*pow(Lmuf2Omur2,2)*QCD::Nc*pow(QCD::Nf,2))/3. + (1460*pow(QCD::Nc,3)*pow(consts::Pi,2))/243. + (4049*consts::egamma*pow(QCD::Nc,3)*pow(consts::Pi,2))/1944. + (469*pow(consts::egamma,2)*pow(QCD::Nc,3)*pow(consts::Pi,2))/108. + (11*pow(consts::egamma,3)*pow(QCD::Nc,3)*pow(consts::Pi,2))/18. + pow(consts::egamma,4)*pow(QCD::Nc,3)*pow(consts::Pi,2) - (1933*Lmuf2Omur2*pow(QCD::Nc,3)*pow(consts::Pi,2))/324. - (77*pow(consts::egamma,2)*Lmuf2Omur2*pow(QCD::Nc,3)*pow(consts::Pi,2))/18. + (121*pow(Lmuf2Omur2,2)*pow(QCD::Nc,3)*pow(consts::Pi,2))/36. + (851*QCD::Nf*pow(consts::Pi,2))/1728. - (Lmuf2Omur2*QCD::Nf*pow(consts::Pi,2))/8. - (26743*pow(QCD::Nc,2)*QCD::Nf*pow(consts::Pi,2))/15552. - (569*consts::egamma*pow(QCD::Nc,2)*QCD::Nf*pow(consts::Pi,2))/1944. - (35*pow(consts::egamma,2)*pow(QCD::Nc,2)*QCD::Nf*pow(consts::Pi,2))/54. - (pow(consts::egamma,3)*pow(QCD::Nc,2)*QCD::Nf*pow(consts::Pi,2))/9. + (1327*Lmuf2Omur2*pow(QCD::Nc,2)*QCD::Nf*pow(consts::Pi,2))/648. + (7*pow(consts::egamma,2)*Lmuf2Omur2*pow(QCD::Nc,2)*QCD::Nf*pow(consts::Pi,2))/9. - (11*pow(Lmuf2Omur2,2)*pow(QCD::Nc,2)*QCD::Nf*pow(consts::Pi,2))/9. + (23*pow(QCD::Nf,2)*pow(consts::Pi,2))/(864.*QCD::Nc) - (539*QCD::Nc*pow(QCD::Nf,2)*pow(consts::Pi,2))/7776. - (10*Lmuf2Omur2*QCD::Nc*pow(QCD::Nf,2)*pow(consts::Pi,2))/81. + (pow(Lmuf2Omur2,2)*QCD::Nc*pow(QCD::Nf,2)*pow(consts::Pi,2))/9. + (29923*pow(QCD::Nc,3)*pow(consts::Pi,4))/77760. - (11*consts::egamma*pow(QCD::Nc,3)*pow(consts::Pi,4))/720. + (43*pow(consts::egamma,2)*pow(QCD::Nc,3)*pow(consts::Pi,4))/180. - (253*Lmuf2Omur2*pow(QCD::Nc,3)*pow(consts::Pi,4))/432. - (11*QCD::Nf*pow(consts::Pi,4))/12960. - (consts::egamma*QCD::Nf*pow(consts::Pi,4))/360. - (277*pow(QCD::Nc,2)*QCD::Nf*pow(consts::Pi,4))/19440. + (23*Lmuf2Omur2*pow(QCD::Nc,2)*QCD::Nf*pow(consts::Pi,4))/216. + (pow(QCD::Nf,2)*pow(consts::Pi,4))/(6480.*QCD::Nc) - (43*QCD::Nc*pow(QCD::Nf,2)*pow(consts::Pi,4))/6480. + (121*pow(QCD::Nc,3)*pow(consts::Pi,6))/6480. - (32707*pow(QCD::Nc,3)*consts::z3)/1296. - (1541*consts::egamma*pow(QCD::Nc,3)*consts::z3)/108. - (88*pow(consts::egamma,2)*pow(QCD::Nc,3)*consts::z3)/9. - 7*pow(consts::egamma,3)*pow(QCD::Nc,3)*consts::z3 + (847*Lmuf2Omur2*pow(QCD::Nc,3)*consts::z3)/108. + (77*consts::egamma*Lmuf2Omur2*pow(QCD::Nc,3)*consts::z3)/6. - (37*QCD::Nf*consts::z3)/8. - (19*consts::egamma*QCD::Nf*consts::z3)/36. - (3*pow(consts::egamma,2)*QCD::Nf*consts::z3)/2. + (11*Lmuf2Omur2*QCD::Nf*consts::z3)/6. + (37*QCD::Nf*consts::z3)/(48.*pow(QCD::Nc,2)) + (5069*pow(QCD::Nc,2)*QCD::Nf*consts::z3)/1296. + (85*consts::egamma*pow(QCD::Nc,2)*QCD::Nf*consts::z3)/54. + (5*pow(consts::egamma,2)*pow(QCD::Nc,2)*QCD::Nf*consts::z3)/18. - (55*Lmuf2Omur2*pow(QCD::Nc,2)*QCD::Nf*consts::z3)/54. - (7*consts::egamma*Lmuf2Omur2*pow(QCD::Nc,2)*QCD::Nf*consts::z3)/3. + (7*pow(QCD::Nf,2)*consts::z3)/(12.*QCD::Nc) - (Lmuf2Omur2*pow(QCD::Nf,2)*consts::z3)/(3.*QCD::Nc) + (5*QCD::Nc*pow(QCD::Nf,2)*consts::z3)/81. + (consts::egamma*QCD::Nc*pow(QCD::Nf,2)*consts::z3)/9. - (2*Lmuf2Omur2*QCD::Nc*pow(QCD::Nf,2)*consts::z3)/27. - (253*pow(QCD::Nc,3)*pow(consts::Pi,2)*consts::z3)/144. - (73*consts::egamma*pow(QCD::Nc,3)*pow(consts::Pi,2)*consts::z3)/36. - (QCD::Nf*pow(consts::Pi,2)*consts::z3)/2. - (13*pow(QCD::Nc,2)*QCD::Nf*pow(consts::Pi,2)*consts::z3)/72. + (3*pow(QCD::Nc,3)*pow(consts::z3,2))/2. + (869*pow(QCD::Nc,3)*consts::z5)/144. + 6*consts::egamma*pow(QCD::Nc,3)*consts::z5 + (5*QCD::Nf*consts::z5)/4. - (5*QCD::Nf*consts::z5)/(4.*pow(QCD::Nc,2)) + (101*pow(QCD::Nc,2)*QCD::Nf*consts::z5)/72.;


    return pow(Lmh2Omuf2,3)*l3c + pow(Lmh2Omuf2,2)*l2c + pow(Lmh2Omuf2,1)*l1c + pow(Lmh2Omuf2,0)*l0c;
  }
  // the different orders of wilson coefficient * xsection expanded in alphas
  double r_wser_0(void){
    return m_LO_();
  }
  double r_wser_1(double Lmh2Omuf2, double Lmuf2Omur2, double w1){
    return m_NLO_(Lmh2Omuf2, Lmuf2Omur2) + 2.0*w1*m_LO_();
  }
  double r_wser_2(double Lmh2Omuf2, double Lmuf2Omur2, double w1, double w2){
    return m_NNLO_(Lmh2Omuf2, Lmuf2Omur2) + 2.0*w1*m_NLO_(Lmh2Omuf2, Lmuf2Omur2) + (pow(w1,2)+2.0*w2)*m_LO_();
  }
  double r_wser_3(double Lmh2Omuf2, double Lmuf2Omur2, double w1, double w2, double w3){
    return m_NNNLO_(Lmh2Omuf2, Lmuf2Omur2) + 2.0*w1*m_NNLO_(Lmh2Omuf2, Lmuf2Omur2) + (pow(w1,2)+2.0*w2)*m_NLO_(Lmh2Omuf2, Lmuf2Omur2)+2.0*(w3+w1*w2)*m_LO_();
  }
  // and the pi^2 resummation corrections
  namespace PI2{
    double r_wser_pi_0(void){
      return 0.0;
    }
    double r_wser_pi_1(void){
      return -(QCD::Nc*pow(consts::Pi,2))/2.;
    }
    double r_wser_pi_2(double Lmh2Omuf2, double Lmuf2Omur2, double w1){
      return (-67*pow(QCD::Nc,2)*consts::Pi)/72. + (5*QCD::Nc*QCD::Nf*consts::Pi)/36. - (r_wser_1(Lmh2Omuf2, Lmuf2Omur2,w1)*QCD::Nc*pow(consts::Pi,2))/2. + (pow(QCD::Nc,2)*pow(consts::Pi,3))/24. + (pow(QCD::Nc,2)*pow(consts::Pi,4))/8.;
    }
    double r_wser_pi_3(double Lmh2Omuf2, double Lmuf2Omur2, double w1, double w2){
      return -(r_wser_2(Lmh2Omuf2, Lmuf2Omur2,w1,w2)*QCD::Nc*pow(consts::Pi,2))/2. + (67*pow(QCD::Nc,3)*pow(consts::Pi,3))/144. - (5*pow(QCD::Nc,2)*QCD::Nf*pow(consts::Pi,3))/72. - (pow(QCD::Nc,3)*pow(consts::Pi,5))/48. - (pow(QCD::Nc,3)*pow(consts::Pi,6))/48. + r_wser_1(Lmh2Omuf2, Lmuf2Omur2,w1)*((-67*pow(QCD::Nc,2)*consts::Pi)/72. + (5*QCD::Nc*QCD::Nf*consts::Pi)/36. + (pow(QCD::Nc,2)*pow(consts::Pi,3))/24. + (pow(QCD::Nc,2)*pow(consts::Pi,4))/8.);
    }
    double r_pi2_exponent(double al){
      return (QCD::Nc*pow(consts::Pi,2)*al)/2. + ((67*pow(QCD::Nc,2)*consts::Pi)/72. - (5*QCD::Nc*QCD::Nf*consts::Pi)/36. - (pow(QCD::Nc,2)*pow(consts::Pi,3))/24.)*pow(al,2);
    }

  }

  // cusp anomalous dimension at different orders
  double r_cusp_A1(void){
    return QCD::CA;
  }
  double r_cusp_A2(void){
    return -(QCD::CA*(-67*QCD::CA + 10*QCD::Nf + 18*QCD::CA*consts::z2))/36.;
  }
  double r_cusp_A3(void){
    return (QCD::CA*(11025*pow(QCD::CA,2) - 2090*QCD::CA*QCD::Nf - 2475*QCD::CF*QCD::Nf - 40*pow(QCD::Nf,2) - 8040*pow(QCD::CA,2)*consts::z2 + 1200*QCD::CA*QCD::Nf*consts::z2 + 2376*pow(QCD::CA,2)*pow(consts::z2,2) + 1980*pow(QCD::CA,2)*consts::z3 - 2520*QCD::CA*QCD::Nf*consts::z3 + 2160*QCD::CF*QCD::Nf*consts::z3))/4320.;
  }
  // 4-loop cusp is unknown, the value here is taken from a pade approximant,
  // so it may be interesting to vary it with a fudge factor
  // numerical dependence of the resummation on A4 is negligible
  double r_cusp_A4(double fudge){
    return (1553*QCD::CA*fudge)/(256.*QCD::CF);
  }
  double r_cusp_D2(void){
    return (QCD::CA*(-202*QCD::CA + 28*QCD::Nf + 198*QCD::CA*consts::z2 - 36*QCD::Nf*consts::z2 + 189*QCD::CA*consts::z3))/54.;
  }
  double r_cusp_D3(void){
    return -(QCD::CA*(297029*pow(QCD::CA,2) - 62626*QCD::CA*QCD::Nf - 46197*QCD::CF*QCD::Nf + 1856*pow(QCD::Nf,2) - 442008*pow(QCD::CA,2)*consts::z2 + 132264*QCD::CA*QCD::Nf*consts::z2 + 11664*QCD::CF*QCD::Nf*consts::z2 - 8640*pow(QCD::Nf,2)*consts::z2 - 541944*pow(QCD::CA,2)*consts::z3 + 100440*QCD::CA*QCD::Nf*consts::z3 + 24624*QCD::CF*QCD::Nf*consts::z3 - 4320*pow(QCD::Nf,2)*consts::z3 + 42768*pow(QCD::CA,2)*consts::z2*consts::z3 + 181764*pow(QCD::CA,2)*consts::z4 - 44712*QCD::CA*QCD::Nf*consts::z4 + 11664*QCD::CF*QCD::Nf*consts::z4 + 139968*pow(QCD::CA,2)*consts::z5))/23328.;
  }
  complex<double> operator*(int n, const complex<double> &x){
    return ((double)n)*x;
  }
  complex<double> operator-(int n, const complex<double> &x){
    return ((double)n)-x;
  }
  complex<double> operator/(int n, const complex<double> &x){
    return ((double)n)/x;
  }
  // the terms in the sudakov exponent
  complex<double> r_g1H(const complex<double> &x){
    return (r_cusp_A1()*(2*x + (1 - 2*x)*log(1 - 2*x)))/(consts::beta_zero*x);
  }
  complex<double> r_g2H(const complex<double> &x, double Lqr, double Lfr){
    return (r_cusp_A2()*(-2*x - log(1 - 2*x)))/pow(consts::beta_zero,2) - (2*consts::egamma*r_cusp_A1()*log(1 - 2*x))/consts::beta_zero + (r_cusp_A1()*(2*Lfr*x + Lqr*log(1 - 2*x)))/consts::beta_zero + (consts::beta_one*r_cusp_A1()*(2*x + log(1 - 2*x) + pow(log(1 - 2*x),2)/2.))/pow(consts::beta_zero,3);
  }
  complex<double> r_g3H(const complex<double> &x, double Lqr, double Lfr){
    return consts::egamma*(-2/consts::beta_zero + 2/(consts::beta_zero*(1 - 2*x)))*r_cusp_A2() + (Lqr/consts::beta_zero - Lqr/(consts::beta_zero*(1 - 2*x)) + (2*Lfr*x)/consts::beta_zero)*r_cusp_A2() + (-1/(2.*pow(consts::beta_zero,2)) + 1/(2.*pow(consts::beta_zero,2)*(1 - 2*x)) - x/pow(consts::beta_zero,2))*r_cusp_A3() + (1/(2.*consts::beta_zero) - 1/(2.*consts::beta_zero*(1 - 2*x)))*r_cusp_D2() + consts::beta_one*r_cusp_A2()*(3/(2.*pow(consts::beta_zero,3)) - 3/(2.*pow(consts::beta_zero,3)*(1 - 2*x)) + x/pow(consts::beta_zero,3) - log(1 - 2*x)/(pow(consts::beta_zero,3)*(1 - 2*x))) + consts::beta_one*consts::egamma*r_cusp_A1()*(2/pow(consts::beta_zero,2) - 2/(pow(consts::beta_zero,2)*(1 - 2*x)) - (2*log(1 - 2*x))/(pow(consts::beta_zero,2)*(1 - 2*x))) + pow(consts::beta_one,2)*r_cusp_A1()*(-1/(2.*pow(consts::beta_zero,4)) + 1/(2.*pow(consts::beta_zero,4)*(1 - 2*x)) - x/pow(consts::beta_zero,4) - log(1 - 2*x)/pow(consts::beta_zero,4) + log(1 - 2*x)/(pow(consts::beta_zero,4)*(1 - 2*x)) + pow(log(1 - 2*x),2)/(2.*pow(consts::beta_zero,4)*(1 - 2*x))) + r_cusp_A1()*(-consts::beta_two/(2.*pow(consts::beta_zero,3)) - 2*pow(consts::egamma,2) - (consts::beta_one*Lqr)/pow(consts::beta_zero,2) + 2*consts::egamma*Lqr - pow(Lqr,2)/2. + consts::beta_two/(2.*pow(consts::beta_zero,3)*(1 - 2*x)) + (2*pow(consts::egamma,2))/(1 - 2*x) + (consts::beta_one*Lqr)/(pow(consts::beta_zero,2)*(1 - 2*x)) - (2*consts::egamma*Lqr)/(1 - 2*x) + pow(Lqr,2)/(2.*(1 - 2*x)) + (consts::beta_two*x)/pow(consts::beta_zero,3) - pow(Lfr,2)*x + (consts::beta_two*log(1 - 2*x))/pow(consts::beta_zero,3) + (consts::beta_one*Lqr*log(1 - 2*x))/(pow(consts::beta_zero,2)*(1 - 2*x)) - 2*consts::z2 + (2*consts::z2)/(1 - 2*x));
  }
  complex<double> r_g4H(const complex<double> &x, double Lqr, double Lfr, double fudge){
    return (-1/(6.*pow(consts::beta_zero,2)) + 1/(6.*pow(consts::beta_zero,2)*pow(1 - 2*x,2)) - (2*x)/(3.*pow(consts::beta_zero,2)))*r_cusp_A4(fudge) + (1/(4.*consts::beta_zero) - 1/(4.*consts::beta_zero*pow(1 - 2*x,2)))*r_cusp_D3() + (2*consts::beta_one*consts::egamma*Lqr*r_cusp_A1()*log(1 - 2*x))/(consts::beta_zero*pow(1 - 2*x,2)) + r_cusp_A3()*((5*consts::beta_one)/(12.*pow(consts::beta_zero,3)) - consts::egamma/consts::beta_zero + Lqr/(2.*consts::beta_zero) - (5*consts::beta_one)/(12.*pow(consts::beta_zero,3)*pow(1 - 2*x,2)) + consts::egamma/(consts::beta_zero*pow(1 - 2*x,2)) - Lqr/(2.*consts::beta_zero*pow(1 - 2*x,2)) + (2*consts::beta_one*x)/(3.*pow(consts::beta_zero,3)) + (2*Lfr*x)/consts::beta_zero - (consts::beta_one*log(1 - 2*x))/(2.*pow(consts::beta_zero,3)*pow(1 - 2*x,2))) + consts::egamma*r_cusp_A2()*(consts::beta_one/pow(consts::beta_zero,2) + 2*Lqr - consts::beta_one/(pow(consts::beta_zero,2)*pow(1 - 2*x,2)) - (2*Lqr)/pow(1 - 2*x,2) - (2*consts::beta_one*log(1 - 2*x))/(pow(consts::beta_zero,2)*pow(1 - 2*x,2))) + r_cusp_D2()*(-consts::beta_one/(4.*pow(consts::beta_zero,2)) + consts::egamma - Lqr/2. + consts::beta_one/(4.*pow(consts::beta_zero,2)*pow(1 - 2*x,2)) - consts::egamma/pow(1 - 2*x,2) + Lqr/(2.*pow(1 - 2*x,2)) + (consts::beta_one*log(1 - 2*x))/(2.*pow(consts::beta_zero,2)*pow(1 - 2*x,2))) + consts::beta_one*r_cusp_A2()*(-Lqr/(2.*pow(consts::beta_zero,2)) + Lqr/(2.*pow(consts::beta_zero,2)*pow(1 - 2*x,2)) + (Lqr*log(1 - 2*x))/(pow(consts::beta_zero,2)*pow(1 - 2*x,2))) + pow(consts::beta_one,2)*r_cusp_A1()*(-(consts::egamma/pow(consts::beta_zero,3)) + Lqr/(2.*pow(consts::beta_zero,3)) - consts::egamma/(pow(consts::beta_zero,3)*pow(1 - 2*x,2)) + Lqr/(2.*pow(consts::beta_zero,3)*pow(1 - 2*x,2)) + (2*consts::egamma)/(pow(consts::beta_zero,3)*(1 - 2*x)) - Lqr/(pow(consts::beta_zero,3)*(1 - 2*x)) + (consts::egamma*pow(log(1 - 2*x),2))/(pow(consts::beta_zero,3)*pow(1 - 2*x,2)) - (Lqr*pow(log(1 - 2*x),2))/(2.*pow(consts::beta_zero,3)*pow(1 - 2*x,2))) + r_cusp_A2()*((-11*pow(consts::beta_one,2))/(12.*pow(consts::beta_zero,4)) + (2*consts::beta_two)/(3.*pow(consts::beta_zero,3)) - 2*pow(consts::egamma,2) - pow(Lqr,2)/2. - pow(consts::beta_one,2)/(12.*pow(consts::beta_zero,4)*pow(1 - 2*x,2)) + consts::beta_two/(3.*pow(consts::beta_zero,3)*pow(1 - 2*x,2)) + (2*pow(consts::egamma,2))/pow(1 - 2*x,2) + pow(Lqr,2)/(2.*pow(1 - 2*x,2)) + pow(consts::beta_one,2)/(pow(consts::beta_zero,4)*(1 - 2*x)) - consts::beta_two/(pow(consts::beta_zero,3)*(1 - 2*x)) - (2*pow(consts::beta_one,2)*x)/(3.*pow(consts::beta_zero,4)) + (2*consts::beta_two*x)/(3.*pow(consts::beta_zero,3)) - 2*pow(Lfr,2)*x + (pow(consts::beta_one,2)*log(1 - 2*x))/(2.*pow(consts::beta_zero,4)*pow(1 - 2*x,2)) + (pow(consts::beta_one,2)*pow(log(1 - 2*x),2))/(2.*pow(consts::beta_zero,4)*pow(1 - 2*x,2)) - 2*consts::z2 + (2*consts::z2)/pow(1 - 2*x,2)) + r_cusp_A1()*((2*pow(consts::beta_one,3))/(3.*pow(consts::beta_zero,5)) - (7*consts::beta_one*consts::beta_two)/(12.*pow(consts::beta_zero,4)) - consts::beta_three/(12.*pow(consts::beta_zero,3)) + (consts::beta_two*consts::egamma)/pow(consts::beta_zero,2) - (4*consts::beta_zero*pow(consts::egamma,3))/3. - (consts::beta_two*Lqr)/(2.*pow(consts::beta_zero,2)) + 2*consts::beta_zero*pow(consts::egamma,2)*Lqr - consts::beta_zero*consts::egamma*pow(Lqr,2) + (consts::beta_zero*pow(Lqr,3))/6. + pow(consts::beta_one,3)/(3.*pow(consts::beta_zero,5)*pow(1 - 2*x,2)) - (5*consts::beta_one*consts::beta_two)/(12.*pow(consts::beta_zero,4)*pow(1 - 2*x,2)) + consts::beta_three/(12.*pow(consts::beta_zero,3)*pow(1 - 2*x,2)) + (consts::beta_two*consts::egamma)/(pow(consts::beta_zero,2)*pow(1 - 2*x,2)) + (4*consts::beta_zero*pow(consts::egamma,3))/(3.*pow(1 - 2*x,2)) - (consts::beta_two*Lqr)/(2.*pow(consts::beta_zero,2)*pow(1 - 2*x,2)) - (2*consts::beta_zero*pow(consts::egamma,2)*Lqr)/pow(1 - 2*x,2) + (consts::beta_zero*consts::egamma*pow(Lqr,2))/pow(1 - 2*x,2) - (consts::beta_zero*pow(Lqr,3))/(6.*pow(1 - 2*x,2)) - pow(consts::beta_one,3)/(pow(consts::beta_zero,5)*(1 - 2*x)) + (consts::beta_one*consts::beta_two)/(pow(consts::beta_zero,4)*(1 - 2*x)) - (2*consts::beta_two*consts::egamma)/(pow(consts::beta_zero,2)*(1 - 2*x)) + (consts::beta_two*Lqr)/(pow(consts::beta_zero,2)*(1 - 2*x)) + (2*pow(consts::beta_one,3)*x)/(3.*pow(consts::beta_zero,5)) - (4*consts::beta_one*consts::beta_two*x)/(3.*pow(consts::beta_zero,4)) + (2*consts::beta_three*x)/(3.*pow(consts::beta_zero,3)) - (consts::beta_one*pow(Lfr,2)*x)/consts::beta_zero + (2*consts::beta_zero*pow(Lfr,3)*x)/3. + (pow(consts::beta_one,3)*log(1 - 2*x))/(2.*pow(consts::beta_zero,5)) - (consts::beta_one*consts::beta_two*log(1 - 2*x))/pow(consts::beta_zero,4) + (consts::beta_three*log(1 - 2*x))/(2.*pow(consts::beta_zero,3)) + (pow(consts::beta_one,3)*log(1 - 2*x))/(2.*pow(consts::beta_zero,5)*pow(1 - 2*x,2)) - (consts::beta_one*consts::beta_two*log(1 - 2*x))/(2.*pow(consts::beta_zero,4)*pow(1 - 2*x,2)) - (2*consts::beta_one*pow(consts::egamma,2)*log(1 - 2*x))/(consts::beta_zero*pow(1 - 2*x,2)) - (consts::beta_one*pow(Lqr,2)*log(1 - 2*x))/(2.*consts::beta_zero*pow(1 - 2*x,2)) - (pow(consts::beta_one,3)*log(1 - 2*x))/(pow(consts::beta_zero,5)*(1 - 2*x)) + (consts::beta_one*consts::beta_two*log(1 - 2*x))/(pow(consts::beta_zero,4)*(1 - 2*x)) - (pow(consts::beta_one,3)*pow(log(1 - 2*x),3))/(6.*pow(consts::beta_zero,5)*pow(1 - 2*x,2)) - 4*consts::beta_zero*consts::egamma*consts::z2 + 2*consts::beta_zero*Lqr*consts::z2 + (4*consts::beta_zero*consts::egamma*consts::z2)/pow(1 - 2*x,2) - (2*consts::beta_zero*Lqr*consts::z2)/pow(1 - 2*x,2) - (2*consts::beta_one*log(1 - 2*x)*consts::z2)/(consts::beta_zero*pow(1 - 2*x,2)) - (8*consts::beta_zero*consts::z3)/3. + (8*consts::beta_zero*consts::z3)/(3.*pow(1 - 2*x,2)));
  }

  // the sudakov exponent
  complex<double> r_sudakov(const complex<double> &n, double as, double Lqr, double Lfr, int order,const string& resummation_type, double fudge  = 1.0){
      complex<double> logN = LogN(n,resummation_type);
    complex<double> rr = 0.0;
    complex<double> lam = as * consts::beta_zero * logN;
    if(order >= 0){
      rr += logN * r_g1H(lam);
    }
    if(order >= 1){
      rr += r_g2H(lam,Lqr,Lfr);
    }
    if(order >= 2){
      rr += as * r_g3H(lam,Lqr,Lfr);
    }
    if(order >= 3){
      rr += pow(as,2) * r_g4H(lam,Lqr,Lfr,fudge);
    }
    return rr;
  }


  complex<double> r_match_LO_(void){
    return -1.0;
  }

  complex<double> r_match_NLO_(const complex<double> &n, double Lqr, double Lfr, unsigned int logOrder, double w1,const string& resummation_type){
      complex<double> logN = LogN(n,resummation_type);
    double a_n0ll = 0.0, a_n1ll = 0.0;//, a_n2ll = 0.0, a_n3ll = 0.0;
    //if(logOrder >= 0)
    a_n0ll = 1.0;
    if(logOrder >= 1)
      a_n1ll = 1.0;
    
    double Lmh2Omuf2 = Lqr -Lfr;
    double Lmuf2Omur2 = Lfr;

    complex<double> res = -r_wser_1(Lmh2Omuf2,Lmuf2Omur2,w1) - 2*r_cusp_A1()*pow(logN,2)*a_n0ll - 2*(2*consts::egamma + Lfr - Lqr)*r_cusp_A1()*logN*a_n1ll;
    return res;
  }
  complex<double> r_match_NNLO_(const complex<double> &n, double Lqr, double Lfr, unsigned int logOrder, double w1, double w2,const string& resummation_type){
      complex<double> logN = LogN(n,resummation_type);
    double a_n0ll = 0.0, a_n1ll = 0.0, a_n2ll = 0.0;
    //if(logOrder >= 0)
    a_n0ll = 1.0;
    if(logOrder >= 1)
      a_n1ll = 1.0;
    if(logOrder >= 2)
      a_n2ll = 1.0;
    
    
    double Lmh2Omuf2 = Lqr -Lfr;
    double Lmuf2Omur2 = Lfr;

    complex<double> res = -r_wser_2(Lmh2Omuf2,Lmuf2Omur2,w1,w2) - 2*pow(r_cusp_A1(),2)*pow(logN,4)*pow(a_n0ll,2) + (-2*r_wser_1(Lmh2Omuf2,Lmuf2Omur2,w1)*(2*consts::egamma + Lfr - Lqr)*r_cusp_A1()*logN - 2*(2*consts::beta_zero*consts::egamma*r_cusp_A1() - consts::beta_zero*Lqr*r_cusp_A1() + r_cusp_A2())*pow(logN,2))*a_n1ll - 2*pow(2*consts::egamma + Lfr - Lqr,2)*pow(r_cusp_A1(),2)*pow(logN,2)*pow(a_n1ll,2) + a_n0ll*(-2*r_wser_1(Lmh2Omuf2,Lmuf2Omur2,w1)*r_cusp_A1()*pow(logN,2) - (4*consts::beta_zero*r_cusp_A1()*pow(logN,3))/3. - 4*(2*consts::egamma + Lfr - Lqr)*pow(r_cusp_A1(),2)*pow(logN,3)*a_n1ll) - logN*(4*consts::beta_zero*pow(consts::egamma,2)*r_cusp_A1() - consts::beta_zero*pow(Lfr,2)*r_cusp_A1() - 4*consts::beta_zero*consts::egamma*Lqr*r_cusp_A1() + consts::beta_zero*pow(Lqr,2)*r_cusp_A1() + 4*consts::egamma*r_cusp_A2() + 2*Lfr*r_cusp_A2() - 2*Lqr*r_cusp_A2() - r_cusp_D2() + 4*consts::beta_zero*r_cusp_A1()*consts::z2)*a_n2ll;
    return res;
  }
  complex<double> r_match_N3LO_(const complex<double> &n, double Lqr, double Lfr, unsigned int logOrder, double w1, double w2, double w3,const string& resummation_type){
    double a_n0ll = 0.0, a_n1ll = 0.0, a_n2ll = 0.0, a_n3ll = 0.0;
    //if(logOrder >= 0)
    a_n0ll = 1.0;
    if(logOrder >= 1)
      a_n1ll = 1.0;
    if(logOrder >= 2)
      a_n2ll = 1.0;
    if(logOrder >= 3)
      a_n3ll = 1.0;
    //double L = -Lqr;
    double Lmh2Omuf2 = Lqr -Lfr;
    double Lmuf2Omur2 = Lfr;

      complex<double> logN = LogN(n,resummation_type);

    complex<double> res = -r_wser_3(Lmh2Omuf2,Lmuf2Omur2,w1,w2,w3) - (4*pow(r_cusp_A1(),3)*pow(logN,6)*pow(a_n0ll,3))/3. + (-2*r_wser_1(Lmh2Omuf2,Lmuf2Omur2,w1)*pow(2*consts::egamma + Lfr - Lqr,2)*pow(r_cusp_A1(),2)*pow(logN,2) - 4*(2*consts::egamma + Lfr - Lqr)*r_cusp_A1()*(2*consts::beta_zero*consts::egamma*r_cusp_A1() - consts::beta_zero*Lqr*r_cusp_A1() + r_cusp_A2())*pow(logN,3))*pow(a_n1ll,2) - (4*pow(2*consts::egamma + Lfr - Lqr,3)*pow(r_cusp_A1(),3)*pow(logN,3)*pow(a_n1ll,3))/3. + pow(a_n0ll,2)*(-2*r_wser_1(Lmh2Omuf2,Lmuf2Omur2,w1)*pow(r_cusp_A1(),2)*pow(logN,4) - (8*consts::beta_zero*pow(r_cusp_A1(),2)*pow(logN,5))/3. - 4*(2*consts::egamma + Lfr - Lqr)*pow(r_cusp_A1(),3)*pow(logN,5)*a_n1ll) + (-(r_wser_1(Lmh2Omuf2,Lmuf2Omur2,w1)*logN*(4*consts::beta_zero*pow(consts::egamma,2)*r_cusp_A1() - consts::beta_zero*pow(Lfr,2)*r_cusp_A1() - 4*consts::beta_zero*consts::egamma*Lqr*r_cusp_A1() + consts::beta_zero*pow(Lqr,2)*r_cusp_A1() + 4*consts::egamma*r_cusp_A2() + 2*Lfr*r_cusp_A2() - 2*Lqr*r_cusp_A2() - r_cusp_D2() + 4*consts::beta_zero*r_cusp_A1()*consts::z2)) - 2*pow(logN,2)*(2*consts::beta_one*consts::egamma*r_cusp_A1() + 4*pow(consts::beta_zero,2)*pow(consts::egamma,2)*r_cusp_A1() - consts::beta_one*Lqr*r_cusp_A1() - 4*pow(consts::beta_zero,2)*consts::egamma*Lqr*r_cusp_A1() + pow(consts::beta_zero,2)*pow(Lqr,2)*r_cusp_A1() + 4*consts::beta_zero*consts::egamma*r_cusp_A2() - 2*consts::beta_zero*Lqr*r_cusp_A2() + r_cusp_A3() - consts::beta_zero*r_cusp_D2() + 4*pow(consts::beta_zero,2)*r_cusp_A1()*consts::z2))*a_n2ll + a_n1ll*(-2*r_wser_2(Lmh2Omuf2,Lmuf2Omur2,w1,w2)*(2*consts::egamma + Lfr - Lqr)*r_cusp_A1()*logN - 2*r_wser_1(Lmh2Omuf2,Lmuf2Omur2,w1)*(2*consts::beta_zero*consts::egamma*r_cusp_A1() - consts::beta_zero*Lqr*r_cusp_A1() + r_cusp_A2())*pow(logN,2) - (4*(consts::beta_one*r_cusp_A1() + 4*pow(consts::beta_zero,2)*consts::egamma*r_cusp_A1() - 2*pow(consts::beta_zero,2)*Lqr*r_cusp_A1() + 2*consts::beta_zero*r_cusp_A2())*pow(logN,3))/3. - 2*(2*consts::egamma + Lfr - Lqr)*r_cusp_A1()*pow(logN,2)*(4*consts::beta_zero*pow(consts::egamma,2)*r_cusp_A1() - consts::beta_zero*pow(Lfr,2)*r_cusp_A1() - 4*consts::beta_zero*consts::egamma*Lqr*r_cusp_A1() + consts::beta_zero*pow(Lqr,2)*r_cusp_A1() + 4*consts::egamma*r_cusp_A2() + 2*Lfr*r_cusp_A2() - 2*Lqr*r_cusp_A2() - r_cusp_D2() + 4*consts::beta_zero*r_cusp_A1()*consts::z2)*a_n2ll) + a_n0ll*(-2*r_wser_2(Lmh2Omuf2,Lmuf2Omur2,w1,w2)*r_cusp_A1()*pow(logN,2) - (4*consts::beta_zero*r_wser_1(Lmh2Omuf2,Lmuf2Omur2,w1)*r_cusp_A1()*pow(logN,3))/3. - (4*pow(consts::beta_zero,2)*r_cusp_A1()*pow(logN,4))/3. + (-4*r_wser_1(Lmh2Omuf2,Lmuf2Omur2,w1)*(2*consts::egamma + Lfr - Lqr)*pow(r_cusp_A1(),2)*pow(logN,3) - (4*r_cusp_A1()*(10*consts::beta_zero*consts::egamma*r_cusp_A1() + 2*consts::beta_zero*Lfr*r_cusp_A1() - 5*consts::beta_zero*Lqr*r_cusp_A1() + 3*r_cusp_A2())*pow(logN,4))/3.)*a_n1ll - 4*pow(2*consts::egamma + Lfr - Lqr,2)*pow(r_cusp_A1(),3)*pow(logN,4)*pow(a_n1ll,2) - 2*r_cusp_A1()*pow(logN,3)*(4*consts::beta_zero*pow(consts::egamma,2)*r_cusp_A1() - consts::beta_zero*pow(Lfr,2)*r_cusp_A1() - 4*consts::beta_zero*consts::egamma*Lqr*r_cusp_A1() + consts::beta_zero*pow(Lqr,2)*r_cusp_A1() + 4*consts::egamma*r_cusp_A2() + 2*Lfr*r_cusp_A2() - 2*Lqr*r_cusp_A2() - r_cusp_D2() + 4*consts::beta_zero*r_cusp_A1()*consts::z2)*a_n2ll) - (logN*(12*consts::beta_one*pow(consts::egamma,2)*r_cusp_A1() + 16*pow(consts::beta_zero,2)*pow(consts::egamma,3)*r_cusp_A1() - 3*consts::beta_one*pow(Lfr,2)*r_cusp_A1() + 2*pow(consts::beta_zero,2)*pow(Lfr,3)*r_cusp_A1() - 12*consts::beta_one*consts::egamma*Lqr*r_cusp_A1() - 24*pow(consts::beta_zero,2)*pow(consts::egamma,2)*Lqr*r_cusp_A1() + 3*consts::beta_one*pow(Lqr,2)*r_cusp_A1() + 12*pow(consts::beta_zero,2)*consts::egamma*pow(Lqr,2)*r_cusp_A1() - 2*pow(consts::beta_zero,2)*pow(Lqr,3)*r_cusp_A1() + 24*consts::beta_zero*pow(consts::egamma,2)*r_cusp_A2() - 6*consts::beta_zero*pow(Lfr,2)*r_cusp_A2() - 24*consts::beta_zero*consts::egamma*Lqr*r_cusp_A2() + 6*consts::beta_zero*pow(Lqr,2)*r_cusp_A2() + 12*consts::egamma*r_cusp_A3() + 6*Lfr*r_cusp_A3() - 6*Lqr*r_cusp_A3() - 12*consts::beta_zero*consts::egamma*r_cusp_D2() + 6*consts::beta_zero*Lqr*r_cusp_D2() - 3*r_cusp_D3() + 12*consts::beta_one*r_cusp_A1()*consts::z2 + 48*pow(consts::beta_zero,2)*consts::egamma*r_cusp_A1()*consts::z2 - 24*pow(consts::beta_zero,2)*Lqr*r_cusp_A1()*consts::z2 + 24*consts::beta_zero*r_cusp_A2()*consts::z2 + 32*pow(consts::beta_zero,2)*r_cusp_A1()*consts::z3)*a_n3ll)/3.;
    return res;
  }
}

ThresholdResummation::ThresholdResummation(const UserInterface &UI, const string &gridDirectory,InputParameters *pInput){
  _input = pInput;
    _resummation_type = UI.giveString("resummation_type");
  Initialize(UI,gridDirectory);
}
ThresholdResummation::~ThresholdResummation(){
  if(_theGrid){
    delete _theGrid;
    _theGrid = NULL;
  }
}

double ThresholdResummation::ResummationCorrection(unsigned int logOrder, unsigned int matchingOrder, bool pi2Resummation){
  double tau = _input->_tau;
  double al = _input->_as_over_pi;
  double ch = Ch(logOrder,false);
  double res = pow(al,2)*_theGrid->DoIntegral([pi2Resummation,tau,logOrder,ch,matchingOrder,this](complex<double> n) -> complex<double>{
    return pow(tau,1.0-n)*(ch*this->SudakovExponential(n,logOrder,pi2Resummation));
  });
  double match = pow(al,2)*_theGrid->DoIntegral([pi2Resummation,tau,logOrder,matchingOrder,this](complex<double> n) -> complex<double>{
    return pow(tau,1.0-n)*(this->MatchingCoefficient(n,matchingOrder,logOrder));
  });
//  cout  << "[ThresholdResummation] logOrder = "<<logOrder
//        << "\tmatching order =  "<<matchingOrder
//        << "\tRes: " << res * _input->_prefactor
//        << "\tmatch: " << match * _input->_prefactor << endl;
  return (match+res) * _input->_prefactor;
}

void ThresholdResummation::Initialize(const UserInterface &UI,
                                      const string &gridDirectory){
    string gname = generate_grid_name(UI);
    if (UI.giveBool("with_resummation_info")) {
        cout << "[ThresholdResummation]: " << "Initializing grid " << gname << endl;
    }
    _theGrid = new GaussGrid(UI,gridDirectory /*+ "/" */+  gname);
    if (UI.giveBool("with_resummation_info")) {
        cout << endl << endl;
        cout << "Ch: " <<  Ch(3,false) << endl;
        cout << "Exp: " << SudakovExponential(complex<double>{2.500000,0.306946},3,false) << endl;
        cout << "Match0: " << MatchingCoefficient(complex<double>{2.500000, 0.306946},0,0) << endl;
        cout << "Match1: " << MatchingCoefficient(complex<double>{2.500000, 0.306946},1,1) << endl;
        cout << "Match2: " << MatchingCoefficient(complex<double>{2.500000, 0.306946},2,2) << endl;
        cout << "Match3: " << MatchingCoefficient(complex<double>{2.500000, 0.306946},3,3) << endl;
    }
    //cout << "[ThresholdResummation]: " << RES::r_sudakov(complex<double>(2.5,1.1), .1,0,0, 3) << endl;
    //cout << "[ThresholdResummation] wc|1: " << RES::r_wser_1(0.0, _input->_wc.c().term_of_order(2).val()) << endl;
    //cout << "[ThresholdResummation] wc|2: " << RES::r_wser_2(0.0, _input->_wc.c().term_of_order(2).val(),_input->_wc.c().term_of_order(3).val()) << endl;
    //cout << "[ThresholdResummation] wc|3: " << RES::r_wser_3(0.0, _input->_wc.c().term_of_order(2).val(),_input->_wc.c().term_of_order(3).val(),_input->_wc.c().term_of_order(4).val()) << endl;
    //cout << "Alpha_s  = " << _input->_as_over_pi * consts::Pi << endl;
    //cout << "Alpha_s/pi  = " << _input->_as_over_pi  << endl;
    //cout << "Testing Ch" << setprecision(8) << endl;
    //double cc = Ch(3,false);
    //cout << "Ch = " << cc << endl;
    //double cp = Ch(3,true);
    //cout << "Ch|pi^2 = " << cp << endl;
    //cout << "exponent = " << RES::PI2::r_pi2_exponent(0.1) << endl;
    //cout << "End testing Ch" << endl;
    //cout << "\tMatching" << endl << endl;
    //cout << "Nlo:\t" << RES::r_match_NLO_(complex<double>(2.5,1.1),0.1,-0.3,3,_input->_wc.c().term_of_order(2).val()) << endl << endl;
    //cout << "Nnlo:\t" << RES::r_match_NNLO_(complex<double>(2.5,1.1),0.1,-0.3,3,_input->_wc.c().term_of_order(2).val(),_input->_wc.c().term_of_order(3).val()) << endl << endl;
    //cout << "N3lo:\t" << RES::r_match_N3LO_(complex<double>(2.5,1.1),0.1,-0.3,3,_input->_wc.c().term_of_order(2).val(),_input->_wc.c().term_of_order(3).val(),_input->_wc.c().term_of_order(4).val()) << endl << endl;
    //cout << "MatchingCoeff:\t" << MatchingCoefficient(complex<double>(2.5,1.3),3,3) << endl;

    //cout << endl << endl;
    //cout << "Expo = NLL" << SudakovExponential(complex<double>(2.5,3.4),1,false) << endl;
    //cout << "Expo = N2LL" << SudakovExponential(complex<double>(2.5,3.4),2,false) << endl;
    //cout << "Expo = N3LL" << SudakovExponential(complex<double>(2.5,3.4),3,false) << endl;
    //cout << "Seri = " << MatchingCoefficient(complex<double>(2.5,3.4),3,3) << endl;
}

double ThresholdResummation::Ch(unsigned int ord, bool pi2){
  double res = 0.0;
  double al = _input->_as_over_pi;
  //double L = _input->_log_muf_over_mh_sq;
  double Lmh2Omuf2 = -_input->_log_muf_over_mh_sq;
  double Lmuf2Omur2 = -_input->_log_mur_over_muf_sq;
  //if(ord >= 0){
    res += RES::r_wser_0();
    if(pi2){
      double dd = RES::PI2::r_wser_pi_0();
      res += pow(al,0)*dd;
    }
  //}
  if(ord >= 1){
    res += al*RES::r_wser_1(Lmh2Omuf2, Lmuf2Omur2,_input->_wc.c().term_of_order(2).val());
    if(pi2){
      double dd = RES::PI2::r_wser_pi_1();
      res += pow(al,1)*dd;
    }
  }
  if(ord >= 2){
    res += pow(al,2)*RES::r_wser_2(Lmh2Omuf2, Lmuf2Omur2,_input->_wc.c().term_of_order(2).val(),_input->_wc.c().term_of_order(3).val());
    if(pi2){
      double dd = RES::PI2::r_wser_pi_2(Lmh2Omuf2, Lmuf2Omur2,_input->_wc.c().term_of_order(2).val());
      res += pow(al,2)*dd;
    }
  }
  if(ord >= 3){
    res += pow(al,3)*RES::r_wser_3(Lmh2Omuf2, Lmuf2Omur2,_input->_wc.c().term_of_order(2).val(),_input->_wc.c().term_of_order(3).val(),_input->_wc.c().term_of_order(4).val());
    if(pi2){
      double dd = RES::PI2::r_wser_pi_3(Lmh2Omuf2, Lmuf2Omur2,_input->_wc.c().term_of_order(2).val(),_input->_wc.c().term_of_order(3).val());
      res += pow(al,3)*dd;
    }
  }
  return res;
}
complex<double> ThresholdResummation::MatchingCoefficient(const complex<double> &n, unsigned int pertOrd, unsigned int logOrd){
  complex<double> res(0.0,0.0);
  double al = _input->_as_over_pi;
  double Lqr = -(_input->_log_muf_over_mh_sq+_input->_log_mur_over_muf_sq);
  double Lfr = -_input->_log_mur_over_muf_sq;
  //if(pertOrd >= 0){
    res += RES::r_match_LO_();
  //}
  if(pertOrd >= 1){
    res += pow(al,1) * RES::r_match_NLO_(n, Lqr, Lfr, logOrd, _input->_wc.c().term_of_order(2).val(),_resummation_type);
  }
  if(pertOrd >= 2){
    res += pow(al,2) * RES::r_match_NNLO_(n, Lqr, Lfr, logOrd, _input->_wc.c().term_of_order(2).val(),_input->_wc.c().term_of_order(3).val(),_resummation_type);
  }
  if(pertOrd >= 3){
    res += pow(al,3) * RES::r_match_N3LO_(n, Lqr, Lfr, logOrd, _input->_wc.c().term_of_order(2).val(),_input->_wc.c().term_of_order(3).val(),_input->_wc.c().term_of_order(4).val(),_resummation_type);
  }
  return res;
}
complex<double> ThresholdResummation::SudakovExponential(const complex<double> &n, unsigned int lOrd, bool pi2Resummation){
  complex<double> se = RES::r_sudakov(n, _input->_as_over_pi, -(_input->_log_muf_over_mh_sq + _input->_log_mur_over_muf_sq), -_input->_log_mur_over_muf_sq, lOrd,_resummation_type);
  double pi2 = exp(RES::PI2::r_pi2_exponent(_input->_as_over_pi));
  if(!pi2Resummation)
    pi2 = 1.0;

    //cout<<"[SudakovExponential] se = "<<se<<endl;
  return exp(se) * pi2;
}

string ThresholdResummation::generate_grid_name(const UserInterface &UI){
  stringstream ss;
  ss << UI.giveString("pdf_set") << '_' << UI.giveInt("pdf_member") << "_" << setiosflags(ios::fixed) << setprecision(3) << UI.giveDouble("muf") << ".bglg";
  return ss.str();
}
