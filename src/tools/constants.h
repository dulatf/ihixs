/**
 * \file  constants.h
 * \brief This file contains constant definitions and shorthands
 */

#ifndef CONSTANTS_H
#define CONSTANTS_H
#include "math.h"
#include <vector>
#include <iostream>

namespace consts
{

    const double Pi=3.1415926535897932385;
    const double pi_square=Pi*Pi;
    const double G_fermi=1.16637e-5;
    const double vev=246.0;
    const double convert_GeV_to_pb=0.389379304*1e9;
    const double convert_GeV_to_fb=0.389379304*1e12;

    const double egamma = 0.57721566490153286061;

    const double z2 = 1.644934066848226;
    const double z3 = 1.202056903159594;
    const double z4 = 1.082323233711138;
    const double z5 = 1.036927755143370;
    const double z6 = 1.017343061984448;

    const double nf=5.0;

    const double beta_zero = 11.0/4.0-consts::nf/6.0;
    const double beta_one = 51.0/8.0-19.0/24.0*consts::nf;
    const double beta_two = 2857.0/128.0-(5033.0/1152.0)*consts::nf
                            +(325.0/3456.0)*pow(consts::nf,2.0);
    const double beta_three = 1.0/256.0*( 149753.0/6.0+3564.0*consts::z3
                                         +(-1078361.0/162.0-6508.0/27.0*consts::z3)*consts::nf
                                         +(50065.0/162.0+6472.0/81.0*consts::z3)*pow(consts::nf,2.0)
                                         +(1093.0/729.0)*pow(consts::nf,3.0));
    // convention: gamma_N = gN/4^N where gN is
    //             the gamma constsant in the pdg quark mass review
    //              http://pdg.lbl.gov/2014/reviews/rpp2014-rev-quark-masses.pdf
    
    const double gamma_one=1;
    const double gamma_two = 4.208333333333333 - 0.1388888888888889* consts::nf;
    
    const double gamma_three =   19.51562500000000
                                - 2.284121493373736 * consts::nf
                                - 0.02700617283950617 * pow(consts::nf,2.);
    
    
    const double gamma_four = 98.94341425520289
                            - 19.10746191863540 * consts::nf
                            + 0.2761625514298946 * pow(consts::nf,2.)
                            + 0.005793222354358383 * pow(consts::nf,3.)
                            ;
    
}

/**
 * \namespace QCD
 * \brief     Namespace containing constants that are characteristic of QCD
 */
namespace QCD
{
    const size_t Nc = 3;            ///< N_c is the number of colors
    const double CA = 3.;           ///< C_A is the Casimir of the adjoint representation
    const long double CF = 4./3.;   ///< C_F is the Casimir of the fundamental representation
    const double TF = 0.5;          ///< T_F is the normalization factor of SU(3) matrices

    const size_t Nf = 5;

    enum Flavor {
        d = 1,
        u = 2,
        s = 3,
        c = 4,
        b = 5,
        t = 6,
        g = 21,
        dbar = -1,
        ubar = -2,
        sbar = -3,
        cbar = -4,
        bbar = -5,
        tbar = -6
    };



}

QCD::Flavor antiParticle(const QCD::Flavor& p);
std::ostream& operator<<(std::ostream& myStream, const QCD::Flavor& flavor);

#endif // CONSTANTS_H
